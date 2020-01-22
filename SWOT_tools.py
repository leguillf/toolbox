#!/usr/bin/env python
# coding: utf-8

import matplotlib
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.cm as cm

import glob
import xarray as xr
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
from collections import Counter

import numpy as np
import GriddedData


grav = 9.81

def gradi(I): 
    """
    Calculates the gradient in the x-direction of an image I and gives as output M.
    In order to keep the size of the initial image the last row is left as 0s.
    """
    
    m, n = I.shape
    M = np.zeros([m,n])

    M[0:-1,:] = np.subtract(I[1::,:], I[0:-1,:])
    return M

def gradj(I): 
    """
    Calculates the gradient in the y-direction of an image I and gives as output M.
    In order to keep the size of the initial image the last column is left as 0s.
    """
    
    m, n = I.shape
    M = np.zeros([m,n])
    M[:,0:-1] =  np.subtract(I[:,1::], I[:,0:-1])
    return M

def div(px, py): 
    """
    Calculates the divergence of a vector (px, py) and gives M as ouput
    ; where px and py have the some size and are the gradient in x and y respectively 
    of a variable p.
    The x component of M (Mx) first row is = to the first row of px.
    The x component of M (Mx) last row is = to - the before last row of px. (last one = 0)
    The y component of M (My) first column is = to the first column of py.
    The y component of M (My) last column is = to - the before last column of py. (last one = 0)
    ??#(de sorte que div=-(grad)^*)
    """
    m, n = px.shape
    M = np.zeros([m,n])
    Mx = np.zeros([m,n])
    My = np.zeros([m,n])
 
    Mx[1:m-1, :] = px[1:m-1, :] - px[0:m-2, :]
    Mx[0, :] = px[0, :]
    Mx[m-1, :] = -px[m-2, :]

    My[:, 1:n-1] = py[:, 1:n-1] - py[:, 0:n-2]
    My[:, 0] = py[:,0]
    My[:, n-1] = -py[:, n-2]
     
    M = Mx + My;
    return M

def laplacian(u):
    """
    Calculates the laplacian of u using the divergence and gradient functions and gives 
    as output Ml.
    """
    Ml = div(gradi(u), gradj(u));
    return Ml

def masked(array):
    """
    """
    m1 = np.ma.masked_invalid(array)
    m2= np.ma.masked_equal(m1,0.)
    return m2

def rel_vort_griddedData(lat, lon, data):
    """
    Calculates the relative vorticity (1/s).  Firstly the laplacian is calculated, and then used to calculate 
    relative vorticity.
    Inputs:
    -
    Output:
    -
    """
    mgrd = GriddedData.grid2D(navlat=lat, navlon=lon)
    corio = GriddedData.corio(mgrd)
    lap = mgrd.div(mgrd.grad(data))
    lap = np.ma.masked_greater(masked(lap),10.)
    
    gx_corio, gy_corio = mgrd.grad(corio)
    beta = gy_corio
    gx_data, gy_data = mgrd.grad(data)
    
    minus_factor = ((grav * beta) / corio**2) * gy_data
    
    rel_vort = (grav * lap / corio ) - minus_factor
    
    return rel_vort

def read_data_nemo(time_index):   
    """
    """
    filenc = datadir_nemo + filename_nemo
    xds = xr.open_dataset(filenc,engine='netcdf4',lock=False)
    if time_index == 'all':
        ssh = xds.sossheig
        tim = xds.time_counter
    elif isinstance(time_index, int):
        ssh = xds.sossheig[time_index,:,:]
        tim = xds.time_counter[time_index]
    else:
        print("Problem")
    lat = xds.nav_lat
    lon = xds.nav_lon
    sshm_hi = ssh.to_masked_array()
    lat_hi = lat.to_masked_array()
    lon_hi = lon.to_masked_array()
    
    return sshm_hi, lat_hi, lon_hi, tim

def plot_natl60(var, vmin, vmax, merv, parv, cmap, ax):
        """
        """
        mapn = Basemap(llcrnrlon=blomin, llcrnrlat=blamin, urcrnrlon=blomax, urcrnrlat=blamax
                          , resolution='i', projection='merc', lat_0 = (blamin+blamax)/2
                          , lon_0 = (blomin+blomax)/2, ax=ax, area_thresh=10)

        mapn.drawcoastlines()
        mapn.fillcontinents(color='#ddaa66')
        mapn.drawmeridians(np.arange(-160, 140, 2), labels=merv, size=18)
        mapn.drawparallels(np.arange(0.5, 70.5, 2), labels=parv, size=18)

        pp = mapn.scatter(lonb, latb, s=20, c=var, linewidth='0', vmin=vmin, vmax=vmax
                         , latlon=True, cmap = cmap)
        pp.set_clim([vmin, vmax])
        
        return pp