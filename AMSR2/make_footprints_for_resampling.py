#     From O:\amsr_l2\v05\resampling\show_antenna_patterns.f


#     The values of these AMSR-E parameters come from E:\AMSR\Resampling_weights\PM1_AMSR_MOD.for
#     The values of th Midori2   parameters come from E:\AMSR\Resampling_weights\ADEOSII_AMSR_MOD.for
#     For these antenna parameters, Peter found one set for 89 GHZ, which was the averaged of 89A and 89B
#     This value was put in the IFREQ=8 slot
#     The IFREQ=6 and 7 slots were set to 0
#     For Midori2, the IFREQ=6 and 7 slots contain the 50.3 and  52.8 patterns respectively.

#     Here, I use slot 7 for 89A and slot 8 for 89B, so they have the same value
#     When I do Midori-2, I will need to decide what to do about the 50 GHz channels (maybe put in average value into slot 6)
  
#  converted to python so I can see what's going on

import numpy as np 

orbit = 2
tdr_path = 'L:/access/resampling/AMSR2/tdr/'
tdr_file = f'{tdr_path}r{orbit:05d}.tdr.extra_fov_1_scan_1.h5'

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

def find_AMSR2_source_locations():

    print(tdr_file)

    d = xr.open_dataset(tdr_file)
    lats = d.latitude.values
    lons = d.longitude.values
    azim = d.azimuth.values

    #for this orbit, the scan places the average lat of all 29x243 scans
    #very close to the equator

    scan0 = 1857
    lats = lats[scan0-28:scan0+29:2,::2]
    lons = lons[scan0-28:scan0+29:2,::2]-lons[scan0,242]
    azim = azim[scan0-28:scan0+29:2,::2]

    return lats,lons,azim

def find_AMSR2_target_locations():

    print(tdr_file)

    d = xr.open_dataset(tdr_file)
    lats = d.latitude.values
    lons = d.longitude.values
    azim = d.azimuth.values

    #for this orbit, the scan places the average lat of all 29x243 scans
    #very close to the equator

    scan0 = 1857
    lats = lats[scan0:scan0+2,:]
    lons = lons[scan0:scan0+2,:]-lons[scan0,242]
    azim = azim[scan0:scan0+2,:]

    return lats,lons,azim

def write_footprint_locations(lats,lons,footprint_type,band = 2):
    import os

    if footprint_type == 'target':
        textfile = 'L:/access/resampling/AMSR2/target_gains_v2/locs/find_target_gain_circular_30km.loc'
    elif footprint_type == 'source':
        textfile = f'L:/access/resampling/AMSR2/source_gains_v2/locs/find_source_gain_AMSR2_band_{band:02d}_locs.txt'
    else:
        raise ValueError(f'Invalid Source Type: {footprint_type}')
    
    sz = lats.shape
    with open(textfile,'w') as f:
        for iscan in range(0,sz[0]):
            for ifov in range(0,sz[1]):
                if footprint_type == 'target':
                    f.write(f' {iscan:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f}\n')
                else:
                    f.write(f' {band:04d} {iscan+1:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f}\n')
                   


    

import xarray as xr 
from matplotlib import pyplot as plt
from matplotlib import cm
import math
import os
import sys
sys.path.append("C:/python_packages_cam/local_earth_grid")
from rss_gridding.local_earth_grid import *
from AMSR2_Antenna_Gain import AMSR2_antenna_gain,target_gain


plt.rcParams.update({'font.size': 18})

freq_names = ['7GHz','11GHz','19GHz','24 GHz','37 GHz']

make_plots = False

footprint_type = 'target'
#footprint_type = 'source'

output_path_binary = 'L:/access/resampling/AMSR2/source_gains_v2/'
output_path_netcdf = 'L:/access/resampling/AMSR2/source_gains_v2_nc/'


#This returns locations with the center lat/lon arrays set to zero by translation
if footprint_type == 'source':
    lats,lons,azim = find_AMSR2_source_locations()
    output_path_binary = 'L:/access/resampling/AMSR2/source_gains_v2/'
    output_path_netcdf = 'L:/access/resampling/AMSR2/source_gains_v2_nc/'
elif footprint_type == 'target':
    lats,lons,azim = find_AMSR2_target_locations()
    output_path_binary = 'L:/access/resampling/AMSR2/target_gains_v2/'
    output_path_netcdf = 'L:/access/resampling/AMSR2/target_gains_v2_nc/'   
else:
    raise ValueError(f'footprint type "{footprint_type}" not valid')

#Write out the locations to text file




#1km spacing very close to Frank's 0.01 degree spacing
grid = LocalEarthGrid(center_lon = 0.0,center_lat = 0.0,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 2401,num_y = 1601,
                          name='Single')

if make_plots:
    grid.addlayer(name='Sum')

x = grid.xlocs*1000.0
y = grid.ylocs*1000.0

#These are nominal values for AMSR2 -- only used for source patterns
slant_range = 1114000.
theta = 55.0

for jfreq in [1,2,3,4]:
    write_footprint_locations(lats,lons,footprint_type,band=jfreq)

    continue
    if jfreq <= 1:
        target_size = 50
    else:
        target_size = 30

    sz = lats.shape
    ifov_end = sz[1]
    iscan_end = sz[0]
        
    for ifov in range(0,ifov_end):
        for iscan in range(0,iscan_end):
            #find the center of footprint in x/y space
            xloc,yloc = grid.projection(lons[iscan,ifov],lats[iscan,ifov])

            print(jfreq,ifov,iscan,lons[iscan,ifov],lats[iscan,ifov],xloc,yloc)

            if footprint_type == 'source':
                # calculate delta theta from boresight vector.
                max_distance_to_consider = 200000.0
                delta_theta,r_norm = calc_delta_theta(slant_range,theta,
                                                azim[iscan,ifov],
                                                x-xloc,y-yloc,
                                                max_distance_to_consider)

                # calculate 2-d antenna gain 
                gain = AMSR2_antenna_gain(delta_theta,jfreq)/np.square(r_norm)

            elif footprint_type == 'target':
                delta_km = 0.001*np.sqrt(np.square(x-xloc)+np.square(y-yloc))
                gain = target_gain(delta_km,diameter_in_km = target_size)

            # store it in the local grid instance
            grid.setdata(gain,name='Single')
            if make_plots:
                grid.adddata(gain,name='Sum')
            
            output_path_binary2 = f'{output_path_binary}band_{jfreq:02d}/'
            os.makedirs(output_path_binary2,exist_ok = True)
            source_file_name = f'{output_path_binary2}s{iscan+1:02d}c{ifov+1:03d}.dat'
            grid.write_binary(outfile = source_file_name,name='Single') 

    #plot it
    if make_plots:
        fig0 = plt.figure(figsize = (10,11))
        outfile = f'{output_path_binary2}plots/{footprint_type}_footprint.png'
        os.makedirs(f'{output_path_binary2}plots/',exist_ok = True)
        grid.contourplot(name='Sum',
                        units='Footprint Weight',
                        title=f'{footprint_type} footprints, {freq_names[jfreq]}',
                        fig=fig0,
                        vmin=0.0,vmax=1.0,
                        cmap = cm.viridis,
                        plt_contours = False,
                        scale=False,
                        cbticks = np.array([0.0,0.2,0.4,0.6,0.8,1.0]),
                        outfile = outfile)
        plt.show()
        print()
    print('Normal End')