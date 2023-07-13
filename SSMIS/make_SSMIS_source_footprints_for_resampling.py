import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
from matplotlib import cm
from read_SSMIS_native_location_file import read_SSMIS_native_location_file
from SSMIS_Antenna_Gain import SSMIS_antenna_gain,target_gain

import sys
sys.path.append('M:/python_packages_cam/rss_gridding/src/rss_gridding/')
from local_earth_grid import LocalEarthGrid

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

def write_footprint_locations(lats,lons,azim,ksat,footprint_type,band = 2):
    import os

    if footprint_type == 'target':
        textfile = Path(f'L:/access/resampling/SSMIS/f{ksat:02d}/target_gains/locs/find_target_gain_f{ksat:02d}.loc')
    elif footprint_type == 'source':
        textfile = Path(f'L:/access/resampling/SSMIS/f{ksat:02d}/source_gains/locs/find_source_gain_f{ksat:02d}_band_{band:02d}_locs.txt')
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
                    print(f' {band:04d} {iscan+1:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f} {azim[iscan,ifov]:12.5f}\n')
                    f.write(f' {band:04d} {iscan+1:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f} {azim[iscan,ifov]:12.5f}\n')
                    


def add_extra_locations(lats,lons,extra_fovs,extra_scans):

    num_native_scans = lats.shape[0]
    num_native_fovs = lats.shape[1]

    num_output_scans = (1 + extra_scans)*num_native_scans - extra_scans
    num_output_fovs  = (1 + extra_fovs)*num_native_scans - extra_fovs

    # set up output arrays
    lats_out = np.full((num_output_scans,num_output_fovs),np.nan,dtype = np.float32)
    lons_out = np.full_like(lats_out,np.nan,dtype=np.float32)
    azim_out = np.full_like(lats_out,np.nan,dtype=np.float32)

    # first add extra fovs to the native scans, if any
    # if extra_fovs > 0:
    #     for jscan in np.arange()
    #     for jfov in np.arange(0,num_output_fovs):
    #         if (jfov % (1+extra_fovs)) == 0:
    #             #just transfer the native location




    # for jscan in np.arange(0,num_output_scans):
    #     for jfov in np.range(0,num_output_fovs):


    # if extra_fovs > 0:
    #     delta_fov = 1.0/extra_fovs

    # if extra_scans > 0:
    #     delta_scan = 1.0/extra_scans

    # if extra_scans == 1:
    #     iscan_list = [0,1]
    # elif extra_scans == 2:
    #     iscan_list = [-1,0,1]
    # elif extra_scans == 3:
    #     iscan_list = [-1,0,1,2]
    # elif extra_scans == 0:
    #     iscan_list = [0]
    # else:
    #     raise ValueError(f'extra_scans = {extra_scans} not supported')


    # for iscan in iscan_list:


    
    # if extra_fovs == 1:
    #     ifov_list = [0,1]
    # elif extra_fovs == 2:
    #     ifov_list = [-1,0,1]
    # elif extra_fovs == 3:
    #     ifov_list = [-1,0,1,2]
    # elif extra_fovs == 0:
    #     ifov_list = [0]
    # else:
    #     raise ValueError(f'extra_scans = {extra_scans} not supported')

    
ksat = 18
orbit_num = 10001
scan0 = 909

freq_list = ['19','22','37']

extra_fovs = 1
extra_scans = 1

make_plots = True

footprint_type = 'source'

locs = read_SSMIS_native_location_file(ksat,orbit_num)

# get the subset of lats and lons
for ifreq,freq in enumerate(freq_list):
    jfreq = ifreq+2
    lats = locs[f'lats{freq}']
    lons = locs[f'lats{freq}']
    azim = locs[f'phi{freq}']
    
    if footprint_type == 'source':
        output_path_binary = Path('L:/access/resampling/SSMIS/source_gains/')
        output_path_netcdf = Path('L:/access/resampling/SSMIS/source_gains_nc/')
        output_path_plot = Path('L:/access/resampling/SSMIS/source_gains_plot/')
    elif footprint_type == 'target':
        output_path_binary = Path(f'L:/access/resampling/SSMIS/target_gains/{target_size:02d}km/')
        output_path_netcdf = Path(f'L:/access/resampling/SSMIS/target_gains_nc/{target_size:02d}km/') 
        output_path_plot = Path(f'L:/access/resampling/SSMIS/target_gains_plot/{target_size:02d}km/') 
    else:
        raise ValueError(f'footprint type "{footprint_type}" not valid')

    os.makedirs(output_path_binary,exist_ok=True)
    os.makedirs(output_path_netcdf,exist_ok=True)
    os.makedirs(output_path_plot,exist_ok=True)

    #find the mean values for slant_range and incidence angle
    rng = locs[f'rng{freq}']
    tht = locs[f'tht{freq}']

    slant_range = np.mean(rng[rng > 10000.0])
    theta = np.mean(tht[tht > 10.0])

    lons = locs[f'lons{freq}']
    lats = locs[f'lats{freq}']
    azim = locs[f'phi{freq}']

    num_cells = lons.shape[1]

    lats = lats[scan0-28:scan0+29,:]
    lons = lons[scan0-28:scan0+29,:]-lons[scan0,int(num_cells/2)]
    azim = azim[scan0-28:scan0+29,:]

    write_footprint_locations(lats,lons,azim,ksat,footprint_type,band = jfreq)
   

    if footprint_type == 'source':
        # for the source footprints, we just use the native locations
        pass
    elif footprint_type == 'target':
        #target locations are the native 37 GHz locations
        lons = locs[f'lons37']
        lats = locs[f'lats37'] 
        # add in the extra locations in between
        lats,lons,azim = add_extra_locations(lats,lons,extra_fovs,extra_scans)
    else:
        raise ValueError(f'footprint type "{footprint_type}" not valid')

    write_footprint_locations(lats,lons,azim,ksat,'target',band = jfreq)

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
            print(freq,jfov,jscan,lons[jscan,jfov],lats[jscan,jfov],xloc,yloc)

            if footprint_type == 'source':
                if (jfov == 0):
                    cross_scan_vector = 0.5*np.array([lons[jscan,1] - lons[jscan,0],
                                                      lats[jscan,1] - lats[jscan,0]])
                elif (jfov == num_fovs-1):
                    cross_scan_vector = 0.5*np.array([lons[jscan,num_fovs-1] - lons[jscan,num_fovs-2],
                                                      lats[jscan,num_fovs-1] - lats[jscan,num_fovs-2]])
                else:
                    cross_scan_vector = 0.25*np.array([lons[jscan,jfov+1] - lons[jscan,jfov-1],
                                                       lats[jscan,jfov+1] - lats[jscan,jfov-1]])

                gain_tot = np.zeros((1601,2401))
                for smear_factor in [-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9]:
                    xloc,yloc = grid.projection(lons[jscan,jfov]+smear_factor*cross_scan_vector[0],
                                                lats[jscan,jfov]+smear_factor*cross_scan_vector[1])

                    # calculate delta theta from boresight vector.
                    max_distance_to_consider = 200000.0
                    delta_theta,r_norm = calc_delta_theta(slant_range,theta,
                                                azim[jscan,jfov],
                                                x-xloc,y-yloc,
                                                max_distance_to_consider)

                    # calculate 2-d antenna gain 
                    gain = SSMIS_antenna_gain(delta_theta,jfreq)/np.square(r_norm)
                    gain_tot = gain_tot + gain
                gain = gain_tot/7.0
            elif footprint_type == 'target':
                delta_km = 0.001*np.sqrt(np.square(x-xloc)+np.square(y-yloc))
                gain = target_gain(delta_km,diameter_in_km = target_size)

            # store it in the local grid instance
            grid.setdata(gain,name='Single')
 
            if make_plots:
                if ((jfov%10 == 3) and (jscan%10 == 0)):
                    grid.adddata(gain,name='Sum')
            
            raise ValueError('stop here')
            output_path_binary2 = f'{output_path_binary}band_{jfreq:02d}/'
            os.makedirs(output_path_binary2,exist_ok = True)
            source_file_name = f'{output_path_binary2}s{jscan+1:02d}c{jfov+1:03d}.dat'
            grid.write_binary(outfile = source_file_name,name='Single',save_thres=0.0001) 

#plot it
    if make_plots:
        fig0 = plt.figure(figsize = (10,11))
        outfile = f'{output_path_binary2}plots/{footprint_type}_footprint.png'
        os.makedirs(f'{output_path_binary2}plots/',exist_ok = True)
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
        png_file = output_path_plot / f'over_scan_{freq}.png'
        fig0.savefig(png_file)
        print()
    print('Normal End')



