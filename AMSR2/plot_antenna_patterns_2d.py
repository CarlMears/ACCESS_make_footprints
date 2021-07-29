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

    orbit = 2
    tdr_path = 'L:/access/resampling/AMSR2/tdr/'
    tdr_file = f'{tdr_path}r{orbit:05d}.tdr.extra_fov_1_scan_1.h5'

    print(tdr_file)

    d = xr.open_dataset(tdr_file)
    lats = d.latitude.values
    lons = d.longitude.values
    azim = d.azimuth.values

    #for this orbit, the scan places the average lat of all 29x243 scans
    #very close to the equator

    scan0 = 1857
    lats = lats[scan0-14:scan0+15,::2]
    lons = lons[scan0-14:scan0+15,::2]-lons[scan0,242]
    azim = azim[scan0-14:scan0+15,::2]


    return lats,lons,azim
    

import xarray as xr 
from matplotlib import pyplot as plt
from matplotlib import cm
import math
import os
import sys
sys.path.append("C:/python_packages_cam/local_earth_grid")
from localearthgrid import *

plt.rcParams.update({'font.size': 18})

output_path_binary = 'L:/access/resampling/AMSR2/source_gains_v2/'
output_path_netcdf = 'L:/access/resampling/AMSR2/source_gains_v2_nc/'

from AMSR2_Antenna_Gain import *

#This returns locations with the center lat/lon arrays set to zero by translation
lats,lons,azim = find_AMSR2_source_locations()

#1km spacing very close to Frank's 0.01 degree spacing
grid = LocalEarthGrid(center_lon = 0.0,center_lat = 0.0,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 2401,num_y = 1601,
                          name='Single')

grid.addlayer(name='Sum')

x = grid.xlocs*1000.0
y = grid.ylocs*1000.0

#These are nominal values for AMSR2
slant_range = 1114000.
theta = 55.0
jfreq = 2
for ifov in range(0,243,20):
    for iscan in range(0,29,14):
        #find the center of footprint in x/y space
        xloc,yloc = grid.projection(lons[iscan,ifov],lats[iscan,ifov])

        print(ifov,iscan,lons[iscan,ifov],lats[iscan,ifov],xloc,yloc)
        # calculate delta theta from boresight vector.
        max_distance_to_consider = 200000.0
        delta_theta,r_norm = calc_delta_theta(slant_range,theta,
                                              azim[iscan,ifov],
                                              x-xloc,y-yloc,
                                              max_distance_to_consider)

        # calculate 2-d antenna gain 
        gain = AMSR2_antenna_gain(delta_theta,jfreq)/np.square(r_norm)

        # store it in the local grid instance
        grid.setdata(gain,name='Single')
        
        output_path_binary2 = f'{output_path_binary}band_{jfreq:02d}/'
        os.makedirs(output_path_binary2,exist_ok = True)
        source_file_name = f'{output_path_binary2}s{iscan+1:02d}c{ifov+1:03d}.dat'
        grid.write_binary(outfile = source_file_name,name='Single') 



#plot it
fig0 = plt.figure(figsize = (10,11))
grid.contourplot(name='AMSR2 Footprints',
                 units='Footprint Weight',
                 fig=fig0,
                 vmin=0.0,vmax=1.0,
                 cmap = cm.viridis,
                 plt_contours = False,
                 scale=False,
                 cbticks = np.array([0.0,0.2,0.4,0.6,0.8,1.0]))
plt.show()
print()