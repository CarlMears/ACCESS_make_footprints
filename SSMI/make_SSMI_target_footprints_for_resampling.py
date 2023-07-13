import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
from matplotlib import cm
from math import floor
from read_SSMI_native_location_file import read_SSMI_native_location_file
from SSMI_Antenna_Gain import SSMI_antenna_gain,target_gain

import sys
sys.path.append('M:/python_packages_cam/rss_gridding/src/rss_gridding/')
from local_earth_grid import LocalEarthGrid

def wt_mean_angle(deg, wts):
    r = np.radians(deg)
    x = np.matmul(wts,np.cos(r))
    y = np.matmul(wts,np.sin(r))
    return np.degrees(np.arctan2(y, x))

def calc_delta_theta_y(range,incidence,dy):
    import numpy as np 
    inc_rad = np.deg2rad(incidence)
    r1 = np.sqrt(range*range + 
                 dy*dy + 
                 2*range*dy*np.sin(np.deg2rad(incidence)))

    z = (np.square(dy) - np.square(r1-range))/(4*range*r1)
    return 2.0*np.rad2deg(np.arcsin(np.sqrt(z))),r1

def calc_delta_theta_x(slant_range,dx):
    import numpy as np
    return np.rad2deg(np.arcsin(dx/slant_range))

def calc_delta_theta(slant_range,theta,azim,x,y,max_distance_to_consider):

    delta_theta = np.ones_like(x)*90.0
    x1 =  np.cos(np.deg2rad(azim))*x - np.sin(np.deg2rad(azim))*y
    y1 =  np.sin(np.deg2rad(azim))*x + np.cos(np.deg2rad(azim))*y

    delta_theta_y,r1 = calc_delta_theta_y(slant_range,theta,y1)
    delta_theta_x = np.rad2deg(np.arcsin(x1/slant_range))

    ok = np.sqrt(np.square(x1)+np.square(y1)) < max_distance_to_consider
    delta_theta[ok] = (np.sqrt(np.square(delta_theta_y) + np.square(delta_theta_x)))[ok]
    return delta_theta,r1/slant_range

def write_footprint_locations(lats,lons,azim,ksat,footprint_type,freq = 2):
    import os

    if footprint_type == 'target':
        textfile = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/target_gains/locs/find_target_gain_f{ksat:02d}_circular.loc')
    elif footprint_type == 'source':
        textfile = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/source_gains/locs/find_source_gain_f{ksat:02d}_freq_{freq:02d}_locs.txt')
    else:
        raise ValueError(f'Invalid Source Type: {footprint_type}')
    
    sz = lats.shape
    os.makedirs(textfile.parent,exist_ok=True)

    with open(textfile,'w') as f:
        for iscan in range(0,sz[0]):
            for ifov in range(0,sz[1]):
                if footprint_type == 'target':
                    f.write(f' {iscan:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f}\n')
                else:
                    print(f' {freq:04d} {iscan+1:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f} {azim[iscan,ifov]:12.5f}\n')
                    f.write(f' {freq:04d} {iscan+1:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f} {azim[iscan,ifov]:12.5f}\n')

def estimate_scan_interp_error(z):

    # estimate the error introduced by interpolating between values of z by
    # calculating the error that occurs if the values at z[n-1],z[n+1] is interpolated
    # to z[n]

    z_interp = 0.5*(z[0:-2] + z[2:])
    z_err_temp = z[1:-1] - z_interp 
    z_err = np.zeros_like(z)
    z_err[1:-1] = z_err_temp
    z_err[0] = z_err[1] - (z_err[2]-z_err[1])
    z_err[-1] = z_err[-2] - (z_err[-3]-z_err[-2])
    return z_err             

def add_extra_locations(lats,lons,azim,extra_fovs,extra_scans):

    num_native_scans = lats.shape[0]
    num_native_fovs = lats.shape[1]

    num_output_scans = (1 + extra_scans)*num_native_scans - extra_scans
    num_output_fovs  = (1 + extra_fovs)*num_native_fovs - extra_fovs

    # set up output arrays
    lats_out = np.full((num_output_scans,num_output_fovs),np.nan,dtype = np.float32)
    lons_out = np.full_like(lats_out,np.nan,dtype=np.float32)
    azim_out = np.full_like(lats_out,np.nan,dtype=np.float32)

    # first add extra fovs to the native scans, if any
    if extra_fovs > 0:
        for jscan in range(0,num_output_scans):
            if jscan % (1+extra_scans) == 0:
                print(f'Adding extra fovs to scan {jscan}')
                jscan_in = int(jscan/(1+extra_scans))
                lon_delta = estimate_scan_interp_error(lons[jscan_in,:])
                lat_delta = estimate_scan_interp_error(lats[jscan_in,:])
                #azim_delta = estimate_scan_interp_error(azim[jscan_in,:]) # too small to bother with
                for jfov in np.arange(0,num_output_fovs):
                    if (jfov % (1+extra_fovs)) == 0:
                        #just transfer the native location
                        jfov_in = int(jfov/(1+extra_fovs))
                        lats_out[jscan,jfov] = lats[jscan_in,jfov_in]
                        lons_out[jscan,jfov] = lons[jscan_in,jfov_in]
                        azim_out[jscan,jfov] = azim[jscan_in,jfov_in]
                    else:
                        jfov_in = int(floor(jfov/(1+extra_fovs)))
                        wt_upper = (jfov - (jfov_in*(1+extra_fovs)))/(1+extra_fovs)
                        wt_lower = 1.0-(wt_upper)
                        error_factor = wt_upper*wt_lower
                        
                        lats_out[jscan,jfov] = wt_lower*lats[jscan_in,jfov_in] + \
                                               wt_upper*lats[jscan_in,jfov_in+1] + \
                                               lat_delta[jscan_in]*error_factor

                        lons_out[jscan,jfov] = wt_lower*lons[jscan_in,jfov_in] + \
                                               wt_upper*lons[jscan_in,jfov_in+1] + \
                                               lon_delta[jscan_in]*error_factor

                        azim_out[jscan,jfov] = wt_mean_angle(azim[jscan_in,jfov_in:jfov_in+2],np.array((wt_lower,wt_upper)))
                        
    #now fill in extra scans
    if extra_scans > 0:
        for jscan in range(1,num_output_scans):
            if jscan % (1+extra_scans) == 0:
                continue
            
            jscan_in = int(jscan/(1+extra_scans))

            wt_upper = (jscan - (jscan_in*(1+extra_scans)))/(1+extra_scans)
            wt_lower = 1.0-(wt_upper)

            jscan_lower = (1+extra_scans)*int(floor(jscan/(1+extra_scans)))
            jscan_upper = jscan_lower + (1+extra_scans)

            print(f'Adding scan {jscan} to output, wt_upper = {wt_upper}, wt_lower = {wt_lower}')
            print(f'combining scans {jscan_lower} and {jscan_upper}')

            lons_out[jscan,:] = wt_mean_angle([lons_out[jscan_lower,:],lons_out[jscan_upper,:]],np.array((wt_lower,wt_upper)))
            lats_out[jscan,:] = wt_mean_angle([lats_out[jscan_lower,:],lats_out[jscan_upper,:]],np.array((wt_lower,wt_upper)))
            azim_out[jscan,:] = wt_mean_angle([azim_out[jscan_lower,:],azim_out[jscan_upper,:]],np.array((wt_lower,wt_upper)))

    azim_out = (azim_out + 1080.0) % 360.0
    return lats_out,lons_out,azim_out


    
ksat = 13
orbit_num = 5000
scan0 = 456
scan0_hi = 912

freq_list = ['19','22','37']

extra_fovs = 3  # add 3 fovs between each native fov
extra_scans = 3 # add 3 scans between each native scan

make_plots = False

footprint_type = 'target'
target_size = 70 #km

locs = read_SSMI_native_location_file(ksat,orbit_num)

# get the subset of lats and lons
for ifreq,freq in enumerate(freq_list):
    
    #all target locations are based on 19 Ghz locations
    lats = locs[f'lats19']
    lons = locs[f'lons19']
    azim = locs[f'phi19']

    delta_scans = 14
    
    output_path_binary = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/target_gains/{target_size:02d}km')
    output_path_plot = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/target_gains/plots') 

    os.makedirs(output_path_binary,exist_ok=True)
    os.makedirs(output_path_plot,exist_ok=True)

    #find the mean values for slant_range and incidence angle

    num_cells = lons.shape[1]

    lats = lats[scan0-delta_scans:scan0+delta_scans+1,:]
    lons = lons[scan0-delta_scans:scan0+delta_scans+1,:]-lons[scan0,int(num_cells/2)]
    azim = azim[scan0-delta_scans:scan0+delta_scans+1,:]


    # add in the extra locations in between the native locations
    lats,lons,azim = add_extra_locations(lats,lons,azim,extra_fovs,extra_scans)
    # the center scan is 56
    # this is one before, and 2 after the center scan,
    # so all three extra scans are included
    lats = lats[55:59,:]
    lons = lons[55:59,:]
    azim = azim[55:59,:]
    write_footprint_locations(lats,lons,azim,ksat,'target',freq = ifreq)
   
    #1km spacing very close to Frank's 0.01 degree spacing
    grid = LocalEarthGrid(center_lon = 0.0,center_lat = 0.0,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 2401,num_y = 1601,
                          name='Single')

    if make_plots:
        grid.addlayer(name='Sum')

    x = grid.xlocs*1000.0
    y = grid.ylocs*1000.0

    num_scans = lats.shape[0]
    num_fovs = lats.shape[1]

    for jscan in np.arange(0,num_scans):
        for jfov in np.arange(0,num_fovs):
        #for jfov in np.arange(0,34):
            xloc,yloc = grid.projection(lons[jscan,jfov],lats[jscan,jfov])
            print(freq,jscan-delta_scans,jfov,lons[jscan,jfov],lats[jscan,jfov],xloc,yloc)

            delta_km = 0.001*np.sqrt(np.square(x-xloc)+np.square(y-yloc))
            gain = target_gain(delta_km,diameter_in_km = target_size)

            # store it in the local grid instance
            grid.setdata(gain,name='Single')
 
            if make_plots:
                print(f'Adding scan {jscan} fov {jfov} to plot')
                grid.adddata(gain,name='Sum')
            
            output_path_binary2 = output_path_binary / f'freq_{ifreq:02d}/'
            os.makedirs(output_path_binary2,exist_ok = True)
            source_file_name = output_path_binary2 /f's{jscan+1:02d}c{jfov+1:03d}.dat'
            grid.write_binary(outfile = source_file_name,name='Single',save_thres=0.0001) 

#plot it
    if make_plots:
        fig0 = plt.figure(figsize = (10,11))
        outfile = output_path_plot / f'decimated_scan_{freq}.png'
        os.makedirs(outfile.parent,exist_ok = True)
        fig0,ax,im =grid.contourplot(name='Sum',
                        units='Footprint Weight',
                        title=f'{footprint_type} footprints, {freq_list[ifreq]}',
                        fig_in=fig0,
                        vmin=0.0,vmax=1.0,
                        cmap = cm.viridis,
                        plt_contours = False,
                        scale=False,
                        cbticks = np.array([0.0,0.2,0.4,0.6,0.8,1.0]),
                        outfile = outfile)
        plt.show()

    print('Normal End')



