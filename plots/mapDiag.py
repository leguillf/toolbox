# Packages

import os
import fnmatch
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
import calendar
import xarray as xr
import matplotlib.pyplot as plt



###############################################################################
# Loading functions
###############################################################################
def load_NATL60(path,file,name_time,name_lon,name_lat,name_var,dt_start,dt_end,dt_ref=datetime(1958,1,1,0,0,0),datetime_type=True):
    """
    NAME 
        load_natl60ssh

    DESCRIPTION
    

        Args:     
        for NATL60 : dt_ref=datetime(1958,1,1,0,0,0)
            
        Param:
        

        Returns: 
        

    """
    # Compute time boundaries in seconds since dt_ref to be compared to timestamps
    time_sec_min = (dt_start - dt_ref).total_seconds()
    time_sec_max = (dt_end - dt_ref).total_seconds() 
    # Read timestamp and grid
    ncin = nc.Dataset(path + file)
    timestamp = np.array(ncin.variables[name_time][:])  
    lon = np.array(ncin.variables[name_lon][:]) 
    lat = np.array(ncin.variables[name_lat][:]) 
    # Find time indexes corresponding to the time boundaries
    idx_time = (timestamp >= time_sec_min) & (timestamp <= time_sec_max)
    out_time = ncin.variables[name_time][idx_time]
    # Read variable 
    var = np.array(ncin.variables[name_var][idx_time,:,:]) 
    

    # Compute datetimes corresponding to the indexes found
    if datetime_type:
        dt_out_time = []
        for itime in range(len(out_time)):
            dt_out_time.append(dt_ref + timedelta(seconds=int(out_time[itime]))) 
    else: 
        dt_out_time = out_time
        
    return var,dt_out_time,lon,lat


def load_DAprod_from_directory(directory, dt_start, dt_end, dt_timestep, name_lon='lat', name_lat='lon', name_var='degraded_sossheig', prefixe = "", suffixe = ""):
    """
    NAME 
        load_prod_from_directory

    DESCRIPTION

        Args:     
        
            
        Param:
        

        Returns: 
        

    """

    listOfFiles = os.listdir(directory)  
    fields = []
    timestamp = []
    dt_curr = dt_start
    first_iter = True
    while dt_curr <= dt_end :          
        yyyy_curr = str(dt_curr.year)
        mm_curr = str(dt_curr.month).zfill(2)
        dd_curr = str(dt_curr.day).zfill(2)
        HH_curr = str(dt_curr.hour).zfill(2)
        MM_curr = str(dt_curr.minute).zfill(2)
        pattern =  prefixe + '_y' + yyyy_curr + 'm' + mm_curr + 'd' + dd_curr + 'h' + HH_curr + MM_curr + suffixe + "*" + ".nc"
        file = [f for f in listOfFiles if fnmatch.fnmatch(f, pattern)]
        if len(file)>1:
            print("Error: several outputs match the pattern... Please set a prefixe.")
            return
        elif len(file)==0:
            print("No file matching the patter : " + pattern)
            dt_curr += dt_timestep
            continue
        else:
            file = file[0]
            fid_deg = nc.Dataset(directory + file) 
            if first_iter:
                lon = np.array(fid_deg.variables[name_lon][:,:]) 
                lat = np.array(fid_deg.variables[name_lat][:,:])  
                first_iter = False
            field_curr = np.mean(fid_deg.variables[name_var][:,:,:],axis=0)
            if np.any(np.isnan(field_curr)):
                print(file,' : NaN value found. Stop here !')
                break
            else: 
                fields.append(field_curr)
                timestamp.append(dt_curr)
        dt_curr += dt_timestep
        
    return np.asarray(fields),np.asarray(timestamp),lon,lat


def load_DUACS(path,file,name_time,name_lon,name_lat,name_var,dt_start,dt_end,dt_ref=datetime(1950,1,1,0,0,0),datetime_type=True):
    ds = xr.open(path + file)
    time_duacs = np.array(ncin.variables[name_time][:]) 
    if 'swot_en_j1_tpn_g2' in file:
        time_duacs += 22919 - 19358
    lon_duacs = np.array(ncin.variables[name_lon][:]) % 360
    lat_duacs = np.array(ncin.variables[name_lat][:]) 
    lon2d_duacs, lat2d_duacs = np.meshgrid(lon_duacs,lat_duacs)
    ssh2d_duacs = np.ma.array(ncin.variables[name_var][:,:,:]) 
    # Convert datetimes to timestamps and select data in the time range
    timestamps_duacs = np.zeros_like(time_duacs)
    for i,dd in enumerate(time_duacs):
        timestamps_duacs[i] = calendar.timegm((dt_ref_duacs + timedelta(days=dd)).timetuple())
        
        
###############################################################################
# Interpolation functions
###############################################################################
    



###############################################################################
# Diagnostic functions
###############################################################################

def diag_multiple_DAfields_comp(refField,expField,lon2d,lat2d,dt_curr,RMSE=None,itime=0,name_RefFields='True state',name_DAfields=None,prefixe='',var='SSH',name_var=None,ncentred=None,var_range=None,cmap='RdBu_r',save=False,path_save=None,plot_err_from_ref=False, xlabel='time (days)',ylabel='RMSE'):
    
    
    Nexp,NY,NX = expField.shape
    if var_range is None:
        if ncentred is not None:
            min_ = expField[:,ncentred:-(ncentred+1),ncentred:-(ncentred+1)].min()
            max_ = expField[:,ncentred:-(ncentred+1),ncentred:-(ncentred+1)].max()
        else:
            min_ = expField.min()        
            max_ = expField.max()
    else:
        min_ = var_range[0]
        max_ = var_range[1]
        
    # Create 2xNda sub plots
    if RMSE is not None:
        gs = gridspec.GridSpec(2, Nexp+2,width_ratios=[1,]*(Nexp+1) + [0.1])
        fig = plt.figure(figsize=(Nexp*10, 5*Nexp)) 
    else:
        gs = gridspec.GridSpec(1, Nexp+2,width_ratios=[1,]*(Nexp+1) + [0.1])
        fig = plt.figure(figsize=(Nexp*10, int(2.5*Nexp)) )
    
    fig.suptitle(dt_curr.strftime("%Y-%m-%d"),fontsize=20)
    
    # Compute error from reference if this option is activated
    if plot_err_from_ref:
        expField = expField - refField
      
    # Plot Reference
    ax0 = plt.subplot(gs[0, 0])
    ax0.pcolormesh(lon2d, lat2d, refField,cmap=plt.cm.get_cmap(cmap),vmin=min_,vmax=max_)
    ax0.set_title(name_RefFields,fontsize=20)

     # Plot subarea where RMSE is computed 
    if ncentred is not None:            
        ax0.plot([lon2d[ncentred,ncentred],lon2d[ncentred,NX-ncentred-1],lon2d[NY-ncentred-1,NX-ncentred-1],lon2d[NY-ncentred-1,ncentred],lon2d[ncentred,ncentred]],
                [lat2d[ncentred,ncentred],lat2d[ncentred,NX-ncentred-1],lat2d[NY-ncentred-1,NX-ncentred-1],lat2d[NY-ncentred-1,ncentred],lat2d[ncentred,ncentred]],'-k')
    # Plot Da fields
    for iexp in range(Nexp):
        ax = plt.subplot(gs[0, iexp+1])
        im = ax.pcolormesh(lon2d, lat2d, expField[iexp],cmap=plt.cm.get_cmap(cmap),vmin=min_,vmax=max_)
        if name_DAfields is not None:
            ax.set_title(name_DAfields[iexp],fontsize=20)
        else:
            ax.set_title('DA product nb ' +str(iexp),fontsize=20)

         # Plot subarea where RMSE is computed 
        if ncentred is not None:            
            ax.plot([lon2d[ncentred,ncentred],lon2d[ncentred,NX-ncentred-1],lon2d[NY-ncentred-1,NX-ncentred-1],lon2d[NY-ncentred-1,ncentred],lon2d[ncentred,ncentred]],
                    [lat2d[ncentred,ncentred],lat2d[ncentred,NX-ncentred-1],lat2d[NY-ncentred-1,NX-ncentred-1],lat2d[NY-ncentred-1,ncentred],lat2d[ncentred,ncentred]],'-k')    
    
    # Colorbar
    cax = plt.subplot(gs[0, -1])
    cbar = fig.colorbar(im, cax=cax, format='%.0e')
    cbar.ax.set_ylabel(var,fontsize=15)
    
    if RMSE is not None:
        Nt = RMSE.shape[1]
        # Plot RMSE 
        ax1 = plt.subplot(gs[1, :])
        for iexp in range(Nexp):
            rmse_t = np.zeros_like(RMSE[0])
            rmse_t[:itime] = RMSE[iexp,:itime]
            rmse_t[itime:] = np.nan
            ax1.plot(rmse_t,label=name_DAfields[iexp])
        ax1.set_xlabel(xlabel,fontsize=15)
        ax1.set_ylabel(ylabel,fontsize=15)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) 
        ax1.legend(fontsize=15)
        ax1.set_ylim([RMSE.min(),RMSE.max()])
        ax1.set_xlim([0,Nt])
        
    # Save figure
    if save:        
        if name_var is not None:
            string = name_var
        else:
            string = var
        title = "comparison_" + string + '_' + prefixe + '_' + dt_curr.strftime("%Y-%m-%d-%H%M")+'.png'
        plt.savefig(path_save + title)
        plt.axis('off')        
        plt.close()
        return title
    else:
        plt.show()   
        return 
    

    



def diag_multiple_DAfields_comp_main(refFields,expFields,lon2d,lat2d,timein,name_RefFields='True state',name_DAfields=None,prefixe='',var='SSH',name_var=None,ncentred=None,var_range=None,plot_rmse=False,cmap='RdBu_r',save=False,path_save=None,plot_err_from_ref=False,run_parallel=False, xlabel='time (days)',ylabel='RMSE'):
    """
    NAME 
        diag_save_2dfields_comp

    DESCRIPTION
    

        Args:     
        
            
        Param:
        

        Returns: 
        

    """
    
    
    if path_save is not None and not os.path.exists(path_save):
        os.makedirs(path_save)
     
    Nt = timein.size
    Nexp = expFields.shape[0]
    
    if var_range is None:
        var_range = [0,0]
        if ncentred is not None:
            var_range[0] = expFields[:,:,ncentred:-(ncentred+1),ncentred:-(ncentred+1)].min()
            var_range[1] = expFields[:,:,ncentred:-(ncentred+1),ncentred:-(ncentred+1)].max()
        else:
            var_range[0] = expFields.min()        
            var_range[1] = expFields.max()
        
    if plot_rmse:
        RMSE = []
        for i in range(Nexp):
            RMSE.append(diag_computing_rmse_soft(refFields, expFields[i], ncentred=ncentred))
        RMSE = np.asarray(RMSE)
    else:
        RMSE = None
    
    if run_parallel:
        from dask import delayed
        ## Defining temporary function for dask parallelization
        def process(iii): 
            diag_multiple_DAfields_comp(
                path_save=path_save,refField=refFields[iii,:,:],expField=expFields[:,iii,:,:],lon2d=lon2d,lat2d=lat2d,itime=iii,dt_curr=timein[iii],RMSE=RMSE,
                name_RefFields=name_RefFields,name_DAfields=name_DAfields,prefixe=prefixe,var=var,ncentred=ncentred,var_range=var_range,cmap=cmap,save=save,plot_err_from_ref=plot_err_from_ref,
                xlabel=xlabel,ylabel=ylabel)
            return 1
 
        # Run diag in parallel on the timeserie
        output = [] 
        for i in range(Nt):
            outpng=delayed(process)(i) 
            output.append(outpng)   
        total = delayed(print)(output)      
        total.compute() 
        
    else:
        for i in range(Nt):
            out = diag_multiple_DAfields_comp(
                    path_save=path_save,refField=refFields[i,:,:],expField=expFields[:,i,:,:],lon2d=lon2d,lat2d=lat2d,itime=i,dt_curr=timein[i],RMSE=RMSE,
                    name_RefFields=name_RefFields,name_DAfields=name_DAfields,prefixe=prefixe,var=var,name_var=name_var,ncentred=ncentred,var_range=var_range,cmap=cmap,save=save,plot_err_from_ref=plot_err_from_ref,
                    xlabel=xlabel,ylabel=ylabel)
            print(out)
    return