#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:23:11 2019

@author: leguillou
"""

import numpy as np 
from scipy import spatial
from scipy.spatial.distance import cdist




def geo2cart(coords):
    """
    NAME 
        ana_transform_coordinates

    DESCRIPTION 
        Transform coordinates from geodetic to cartesian
        
        Args: 
            coords : a set of lan/lon coordinates (e.g. a tuple or 
             an array of tuples)


        Returns: a set of cartesian coordinates (x,y,z)
            
    """ 

    # WGS 84 reference coordinate system parameters
    A = 6378.137 # major axis [km]   
    E2 = 6.69437999014e-3 # eccentricity squared    

    coords = np.asarray(coords).astype(np.float)

    # is coords a tuple? Convert it to an one-element array of tuples
    if coords.ndim == 1:
        coords = np.array([coords])

    # convert to radiants
    lat_rad = np.radians(coords[:,1])
    lon_rad = np.radians(coords[:,0]) 

    # convert to cartesian coordinates
    r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
    x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
    y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
    z = r_n * (1 - E2) * np.sin(lat_rad)

    return np.column_stack((x, y, z))



def construct_ground_pixel_tree(lon,lat):    
    coords = np.column_stack((lon.ravel(),lat.ravel()))
    # construct KD-tree
    ground_pixel_tree = spatial.cKDTree(geo2cart(coords))
    subdomain = geo2cart(coords)[0:100]
    eucl_dist = cdist(subdomain, subdomain, metric="euclidean")
    dist_threshold = np.min(eucl_dist[np.nonzero(eucl_dist)])
    
    return ground_pixel_tree, dist_threshold
   
    
    

def project_obsvar_to_state_grid(var,lon,lat,ground_pixel_tree,n_neighbours,dist_threshold,ny,nx,dist_scale):
    """
    NAME 
        ana_project_nadir_obs_to_state_grid

    DESCRIPTION 
        Project the observations to the state grid by seeking the nearest grid point 
        
        Args: 
            sat_info (Satellite object): information specific to the satellite currently processed 
            obs_file (strings): path and name of observation
            state_vector (string): path and names of state_vector
          
        Param: from swotda_params_specific.py

        Returns: 
             obs_projected (mask array) : observation projected on the state grid 
             filter_projected (mask array) : coefficient for filtering the 2D data for nudging
                 
            
    """     
     
    ###########
    # Reshaping
    ########## 
    lon = lon.ravel()
    lat = lat.ravel()
    var = var.ravel()
        
    ##########
    # Projecting the observations on the model grid
    ##########
    dict_var = {}
    dict_dist = {}
    if var.size>0:
        for i_obs,(lon_,lat_) in enumerate(zip(lon,lat)):                
            dist,ind = ground_pixel_tree.query(geo2cart((lon_,lat_)),k=n_neighbours)  
            if n_neighbours>1:
                dist = dist[0]
                ind = ind[0]
            for k in range(n_neighbours):
                if np.isfinite(dist[k]) and np.isfinite(var[i_obs]): 
                    ind_ = ind[k]
                    if ind_ in dict_var:
                        dict_var[ind_] = np.append(dict_var[ind_],var[i_obs])
                        dict_dist[ind_] = np.append(dict_dist[ind_],dist[k])
                    else:
                        dict_var[ind_] = np.array([var[i_obs]])
                        dict_dist[ind_] = np.array([dist[k]])

                                                       
    # Average on each pixel
    obs_projected = np.zeros(ny*nx)
    obs_projected[:] = np.nan
    tapering = np.zeros(ny*nx)
    tapering[:] = np.nan
    for pix in dict_var:
        obs_inside_pix = np.where((dict_dist[pix]<dist_threshold))[0]
        if np.any(obs_inside_pix):
            obs_projected[pix] = np.average(dict_var[pix][obs_inside_pix],weights=np.exp(-(dict_dist[pix][obs_inside_pix]**2/(0.5*dist_scale)**2)))
            tapering[pix] = 1
        else:
            obs_projected[pix] = np.average(dict_var[pix],weights=np.exp(-(dict_dist[pix]**2/(0.5*dist_scale)**2)))
            tapering[pix] = np.exp(-np.min(dict_dist[pix]/dist_scale))
    # Normalize between 0 and 1:
    tapering = (tapering-np.nanmin(tapering))/(np.nanmax(tapering)-np.nanmin(tapering))
        
    # Mask the useless pixels
    obs_projected = np.ma.masked_invalid(obs_projected).reshape(ny,nx) 
    obs_projected = np.ma.masked_greater(obs_projected, 50)
    tapering = np.ma.masked_invalid(tapering).reshape(ny,nx)                       
            
    return obs_projected,tapering



                                                       
            
    return obs_projected

