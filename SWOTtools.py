#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:20:07 2019

@author: leguillou
"""



class swot_pass():
    def __init__(self, datadir, namfiles, namdescr):
        self.datadir = datadir
        self.namfiles = namfiles
        self.namdescr = namdescr
        fdescr = open(self.datadir + self.namdescr)
        flines = fdescr.readlines()
        fdescr.close()
        exec('self.' + flines[5]) # corresponds to delta_ac [km]
        exec('self.' + flines[6]) # corresponds to delta_al [km]
        exec('self.' + flines[16]) # corresponds to halfgap size [km] 
        
        self.dx = self.delta_ac * 1000 # so in m
        self.dy = self.delta_al * 1000 # so in m
        self.grid_size = self.delta_ac * 1000 # so in m # as regular SWOT grid, can use either delta_ac or _al
        
    def get_list_of_swotfiles(self):
        """
        Obtain a list of all the available SWOT files
        Function(s) used in:
        - get_list_of_cycles
        - get_list_of_passes
        """
        return glob.glob(self.datadir + self.namfiles + '*.nc')

    def get_list_of_cycles(self):
        """
        Obtain a list of all the possible cycles
        """
        swotfiles = self.get_list_of_swotfiles()
        swotcycles = [s.split('_c')[-1].split('_')[0] for s in swotfiles]
        npcycles = np.array(swotcycles,dtype=int)
        # To get back the cycle number as a list instead of array:    
        dictcount = Counter(npcycles)
        return dictcount.keys()

    def get_list_of_passes(self):
        """
        Obtain a list of all the possible passes
        """
        swotfiles = self.get_list_of_swotfiles()
        swotpasses = [s.split('_p')[-1].split('.')[0] for s in swotfiles]
        nppasses = np.array(swotpasses,dtype=int)
        # To get back the cycle number as a list instead of array:
        dictcount = Counter(nppasses)
        passes = dictcount.keys()
        passes.sort()
        return passes
    
    def read_data_nc(self, ncycle, npass):
        """
        Read SWOT data from netcdf file for a particular:
        -ncycle: Cycle number (run get_list_of_cycles)
        -npass: Pass number (run get_list_of_passes)
        """
        self.myfile = self.datadir + self.namfiles \
                + str(ncycle).zfill(2) \
                + '_p' + str(npass).zfill(3) + '.nc'
        nc = Dataset(self.myfile)
        self.lon = nc.variables['lon'][:]
        self.lat = nc.variables['lat'][:]
        self.time = nc.variables['time'][:]
        self.SSH_obs = nc.variables['SSH_obs'][:]
        self.SSH_model = nc.variables['SSH_model'][:]
        nc.close()    
        
        self.nhalfswath = np.shape(self.lon)[1]/2
        
        
    def plot_single_pass(self, lon, lat, var, vmin, vmax, merv, parv, cmap, **options):
        """
        Plots an individual plot of the data using scatter and Basemap.
        Input(s):
        - var = variable to be plotted
        - box = box region of the area wanted to be shown
        , where box is a 1 x 4 array: 
        [minimum_longitude maximum_longitude minimum_latitude maximum_latitude]
        - vmin = minimum value of the colorbar 
        - vmax = maximum value of the colorbar
        - merv = a 1 x 4 array to know if to label and where the meridians
        -- [0 0 0 0] = no labels
        -- [1 0 0 1]
        - parv = like merv, but for the parallels' labels
        - cmap = colormap to be used
        - options:
            -- ax = axis (for e.g. if only one plot, not necessary)
            -- extend_opt = 'both', 'max', or 'min', to extend colorbar
        Output(s):
        - pp = plot object
        """
        lomin = swot_pass1.box[0]
        lomax = swot_pass1.box[1]
        lamin = swot_pass1.box[2]
        lamax = swot_pass1.box[3]

        if not options:
            map = Basemap(llcrnrlon=lomin, llcrnrlat=lamin, urcrnrlon=lomax, urcrnrlat=lamax
                          , resolution='i', projection='merc', lat_0 = (lamin+lamax)/2
                          , lon_0 = (lomin+lomax)/2, area_thresh=10)
        else:
            map = Basemap(llcrnrlon=lomin, llcrnrlat=lamin, urcrnrlon=lomax, urcrnrlat=lamax
                          , resolution='i', projection='merc', lat_0 = (lamin+lamax)/2
                          , lon_0 = (lomin+lomax)/2, ax=options.get("ax"), area_thresh=10)

        map.drawcoastlines()

        map.fillcontinents(color='#ddaa66')
        map.drawmeridians(np.arange(-160, 140, 2), labels=merv, size=18)
        map.drawparallels(np.arange(0.5, 70.5, 2), labels=parv, size=18)

        pp = map.scatter(lon, lat, s=2, c=var, linewidth='0', vmin=vmin, vmax=vmax, latlon=True
                         , cmap = cmap)
        pp.set_clim([vmin, vmax])

        if not options:
            cbar = plt.colorbar(pp, extend=extend_opt)
            cbar.ax.tick_params(labelsize=16)
            
        return pp
    
    
    
    