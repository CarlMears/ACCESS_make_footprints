import numpy as np 
import xarray as xr 
from matplotlib import pyplot as plt
from matplotlib import cm
from pathlib import Path
from read_SSMI_native_location_file import read_SSMI_native_location_file

import math
import sys
sys.path.append('M:/python_packages_cam/rss_gridding/src/rss_gridding/')
from local_earth_grid import LocalEarthGrid

plt.rcParams.update({'font.size': 18})
#show footprint Locations

wt_path = Path('L:/access/resampling/SSMI/F13/resample_weights/70km')
freq = 19
ifreq = 0
scan = 1
center_fov = 243
if freq == 19:
    major = 69
    minor = 43
elif freq == 22:
    major = 50
    minor = 40
elif freq == 37:
    major = 37
    minor = 28
elif freq == 85:
    major = 15
    minor = 13
else:
    raise ValueError(f'Freq: {freq} not valid')

num_source_fovs = 64
num_source_scans = 29
num_target_fovs = 253
num_target_scans = 4 

wt_file = wt_path / f'resample_weights_band_{ifreq:02d}.dat'
temp = np.fromfile(wt_file,dtype=np.float64)
wts = np.reshape(temp,(num_target_scans,num_target_fovs,num_source_scans,num_source_fovs))

mx_wts = np.nanmax(wts)
scale_max = (1+np.floor(10.0*mx_wts))*0.1
print(mx_wts,scale_max)

locs = read_SSMI_native_location_file(13,5000)
freq = '19'
lats = locs[f'lats{freq}']
lons = locs[f'lons{freq}']
azim = locs[f'phi{freq}']

delta_scans = 14
scan0 = 456

num_cells = lons.shape[1]

lats = lats[scan0-delta_scans:scan0+delta_scans+1,:]
lons = lons[scan0-delta_scans:scan0+delta_scans+1,:]-lons[scan0,int(num_cells/2)]
azim = azim[scan0-delta_scans:scan0+delta_scans+1,:]



center_lon = 0.0
center_lat = 0.0

grid = LocalEarthGrid(center_lon = 0.0,center_lat = 0.0,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 2401,num_y = 1601,
                          name='SSMI Footprints')
major=5
minor=3
for scan in range(0,29):
    for fov in range(0,64):
        lon = lons[scan,fov]
        lat = lats[scan,fov]
        grid.addgaussian_at_lonlat(lon,lat,
                            major,minor,90.0+azim,0.75*scale_max,normalize = False,
                            name = 'AMSR2 Footprints')

fig0 = plt.figure(figsize = (10,11))
grid.contourplot(name='SSMI Footprints',
                 units='Footprint Weight',
                 fig=fig0,
                 vmin=-scale_max,vmax=scale_max,
                 cmap = cm.seismic,
                 plt_contours = False,
                 scale=False,
                 cbticks = scale_max*np.array([-1.0,-0.5,0.0,0.5,1.0]))

plt.show()
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



