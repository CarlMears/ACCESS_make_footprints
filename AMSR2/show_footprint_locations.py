import numpy as np 
import xarray as xr 
from matplotlib import pyplot as plt
from matplotlib import cm

import math
import sys
sys.path.append("C:/python_packages_cam/local_earth_grid")
from localearthgrid import *

plt.rcParams.update({'font.size': 18})
#show footprint Locations

wt_path = 'L:/access/resampling/AMSR2/resample_weights/circular_30km/netcdf/'
freq = 19
scan = 1
center_fov = 243
if freq == 6:
    major = 62
    minor = 35
elif freq == 7:
    major = 62
    minor = 35
elif freq == 11:
    major = 42
    minor = 24
elif freq == 19:
    major = 22
    minor = 14
elif freq == 24:
    major = 19
    minor = 11
elif freq == 37:
    major = 12
    minor = 7
elif freq == 89:
    major = 5
    minor = 3
else:
    raise ValueError(f'Freq: {freq} not valid')


wt_file = f'{wt_path}resample_weight_slice_scan_{freq:02d}_{scan:03d}_fov_{center_fov:03d}.nc'
d = xr.open_dataset(wt_file)
wts = d.wts.values
mx_wts = np.nanmax(wts)
scale_max = (1+np.floor(10.0*mx_wts))*0.1
print(mx_wts,scale_max)

orbit = 2
tdr_path = 'L:/access/resampling/AMSR2/tdr/'
tdr_file = f'{tdr_path}r{orbit:05d}.tdr.extra_fov_1_scan_1.h5'
print(tdr_file)
d = xr.open_dataset(tdr_file)
lats = d.latitude.values
lons = d.longitude.values
azi = d.azimuth.values

print(lats[1903,242])
print(lons[1903,242])

center_lon = 184.51532
center_lat = -0.0049385205
center_scan = 1813
center_fov = 242

# print(lats[center_scan,center_fov])
# print(lons[center_scan,center_fov])
grid = LocalEarthGrid(center_lon = center_lon,center_lat = center_lat,
                          delta_x = 0.2,delta_y = 0.2,
                          num_x = 401,num_y = 401,
                          name='AMSR2 Footprints')

 grid.addgaussian_at_lonlat(center_lon,center_lat,
                            30.0,30.0,0.0,0.75*scale_max,normalize = False,
                            name = 'AMSR2 Footprints')

 lat = lats[center_scan-8,center_fov+14]
 lon = lons[center_scan-8,center_fov+14]
 azim = azi[center_scan-8,center_fov+14]

 grid.addgaussian_at_lonlat(lon,lat,
                            major,minor,90.0+azim,0.75*scale_max,normalize = False,
                            name = 'AMSR2 Footprints')

fig0 = plt.figure(figsize = (10,11))
grid.contourplot(name='AMSR2 Footprints',
                 units='Footprint Weight',
                 fig=fig0,
                 vmin=-scale_max,vmax=scale_max,
                 cmap = cm.seismic,
                 plt_contours = False,
                 scale=False,
                 cbticks = scale_max*np.array([-1.0,-0.5,0.0,0.5,1.0]))

 png_file = f'C:/job_access/resampling/prelim_results/Footprints/Footprint_locs_patterns_only_{freq:02d}.png'
fig0.savefig(png_file)
 grid.plot_ellipse(0.0,0.0,
                 15.,15.,
                 0.0,
                 color='black',alpha = 1.0,fig=fig0,facecolor=plt.cm.BrBG(0.75))

major = 22.0
minor = 14.0
ax = None


for scan in range(center_scan-28,center_scan+29):
    for fov in range(center_fov-50,center_fov+51):
        wt_scan = int((scan-center_scan)/2.0)+14
        wt_fov = int(fov/2)
        print(wt_scan,wt_fov,wts[wt_scan,wt_fov])
        lat = lats[scan,fov]
        lon = lons[scan,fov]
        azim = azi[scan,fov]
        print(lat)
        x,y = grid.projection(lon,lat)
        x = x/1000.0
        y = y/1000.0
            #print overall_grid.x_extent
            #print overall_grid.y_extent
        if ((x > grid.x_extent[0]) and (x < grid.x_extent[1]) and 
            (y > grid.y_extent[0]) and (y < grid.y_extent[1])):
            print(x,y,lon,lat) #,footprint['tb'],old_footprint_list_15[i]['tb']
            scaled_wt = 0.5 + (0.5/scale_max)*wts[wt_scan,wt_fov]
            if scaled_wt < 0.0:
                scaled_wt = 0.0
            if scaled_wt > 1.0:
                 scaled_wt = 1.0
            
            #facecolor = plt.cm.seismic(scaled_wt)
            facecolor='red'
            grid.plot_ellipse(x,y,
                1.4,1.4,
                azim,
                color='black',alpha = 1.0,fig=fig0,facecolor=facecolor,ax=ax)



for scan in range(center_scan-28,center_scan+29,2):
    for fov in range(center_fov-50,center_fov+51,2):
        wt_scan = int((scan-center_scan)/2.0)+14
        wt_fov = int(fov/2)
        print(wt_scan,wt_fov,wts[wt_scan,wt_fov])
        lat = lats[scan,fov]
        lon = lons[scan,fov]
        azim = azi[scan,fov]
        print(lat)
        x,y = grid.projection(lon,lat)
        x = x/1000.0
        y = y/1000.0
            #print overall_grid.x_extent
            #print overall_grid.y_extent
        if ((x > grid.x_extent[0]) and (x < grid.x_extent[1]) and 
            (y > grid.y_extent[0]) and (y < grid.y_extent[1])):
            print(x,y,lon,lat) #,footprint['tb'],old_footprint_list_15[i]['tb']
            scaled_wt = 0.5 + (0.5/scale_max)*wts[wt_scan,wt_fov]
            if scaled_wt < 0.0:
                scaled_wt = 0.0
            if scaled_wt > 1.0:
                 scaled_wt = 1.0
            
            #facecolor = plt.cm.seismic(scaled_wt)
            facecolor='blue'
            grid.plot_ellipse(x,y,
                major/20.0,minor/20.0,
                azim,
                color='black',alpha = 1.0,fig=fig0,facecolor=facecolor,ax=ax)
png_file = f'C:/job_access/resampling/prelim_results/Footprints/Footprint_locs_schematic.png'
fig0.savefig(png_file)
plt.show()
print



