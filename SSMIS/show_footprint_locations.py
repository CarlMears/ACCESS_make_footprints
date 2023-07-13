import sys

sys.path.append("M:/job_access/python/land_water/")
sys.path.append("M:/job_access/python/resample_wts/")

import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from common.read_gain_array import read_gain_array
sys.path.append("C:/python_packages_cam/local_earth_grid")
from rss_gridding.local_earth_grid import LocalEarthGrid
from common.random_land_mask import SomewhatRandomLandWater

import matplotlib.pyplot as plt
from numpy.random import RandomState
from numpy.random import default_rng

from scipy.ndimage import shift
from rss_gridding.quadrilateral_interpolation import InterpBilinearQuadrilateral



def read_target_location_file(sat_name='SSMIS',ksat=18,beamwidth=70):

    path = f"L:access/resampling/{sat_name}/f{ksat:02d}/target_gains/locs/"
    file_name = f"{path}find_target_gain_f{ksat:02d}_circular.loc"

    scan = []
    fov = []
    lats = []
    lons = []
    with open(file_name,'r') as f:
        for line in f:
            x = line.split()
            scan.append(int(x[0]))
            fov.append(int(x[1]))
            lats.append(float(x[2]))
            lons.append(float(x[3]))
           
    scan = np.array(scan)
    fov = np.array(fov)
    lats = np.array(lats)
    lons = np.array(lons)

    num_target_scans = np.max(scan)-np.min(scan)+1
    num_target_fovs = np.max(fov)-np.min(fov) + 1

    print(num_target_scans,num_target_fovs,num_target_scans*num_target_fovs,len(scan))

    assert(num_target_scans*num_target_fovs == len(scan))
    print()

    scan = np.reshape(scan,(num_target_scans,num_target_fovs))
    fov = np.reshape(fov,(num_target_scans,num_target_fovs))
    lats = np.reshape(lats,(num_target_scans,num_target_fovs))
    lons = np.reshape(lons,(num_target_scans,num_target_fovs))

    return lats,lons

def read_source_location_file(sat_name='SSMIS',ksat=18,band_num=4):

    path = f"L:access/resampling/{sat_name}/f{ksat:02d}/source_gains/locs/"
    file_name = f"{path}find_source_gain_f{ksat:02d}_band_{band_num:02d}_locs.txt"

    band=[]
    scan = []
    fov = []
    lats = []
    lons = []
    azim = []
    with open(file_name,'r') as f:
        for line in f:
            x = line.split()
            band.append(int(x[0]))
            scan.append(int(x[1]))
            fov.append(int(x[2]))
            lats.append(float(x[3]))
            lons.append(float(x[4]))
            azim.append(float(x[5]))

    num_source_scans = np.max(scan)-np.min(scan)+1
    num_source_fovs = np.max(fov)-np.min(fov) + 1
    print(num_source_scans,num_source_fovs,num_source_scans*num_source_fovs,len(scan))
    assert(num_source_scans*num_source_fovs == len(scan))   
    scan = np.array(scan)
    fov = np.array(fov)
    lats = np.array(lats)
    lons = np.array(lons)
    azim = np.array(azim)

    # scan = np.reshape(scan,(num_source_scans,num_source_fovs))
    # fov = np.reshape(fov,(num_target_scans,num_target_fovs))
    lats = np.reshape(lats,(num_source_scans,num_source_fovs))
    lons = np.reshape(lons,(num_source_scans,num_source_fovs))
    azim = np.reshape(azim,(num_source_scans,num_source_fovs))

    return lats,lons,azim

center_lon = 0.02
center_lat = 2.06

sensor_name = 'SSMIS'
ksat = 18

plt.rcParams.update({'font.size': 18})

grid = LocalEarthGrid(center_lon = center_lon,center_lat = center_lat,
                          delta_x = 0.03,delta_y = 0.03,
                          num_x = 401,num_y = 401,
                          name='SSMIS Footprints')

target_lats,target_lons = read_target_location_file(sat_name=sensor_name,ksat=ksat,beamwidth=70)
source_lats,source_lons,source_azim = read_source_location_file(sat_name=sensor_name,ksat=ksat,band_num=4)

center_scan = 28
major = 64
minor = 28
fig0 = plt.figure(figsize = (10,11))
ax = None
for scan in [0,1]:
    for fov in range(0,446):
        lat = target_lats[scan,fov]
        lon = target_lons[scan,fov]
        x,y = grid.projection(lon,lat)
        x = x/1000.0
        y = y/1000.0
        if ((x > grid.x_extent[0]) and (x < grid.x_extent[1]) and 
            (y > grid.y_extent[0]) and (y < grid.y_extent[1])):
            print(x,y,lon,lat) #,footprint['tb'],old_footprint_list_15[i]['tb']
            facecolor='red'
            ax = grid.plot_ellipse(x,y,
                major/250.0,major/250.0,                 
                0.0,
                color='black',alpha = 1.0,fig=fig0,facecolor=facecolor,ax=ax)

# for scan in range(-28,29):
#     for fov in range(0,90):
#         lat = source_lats[scan+center_scan,fov]
#         lon = source_lons[scan+center_scan,fov]
#         azim = source_azim[scan+center_scan,fov]
#         x,y = grid.projection(lon,lat)
#         x = x/1000.0
#         y = y/1000.0
#         if ((x > grid.x_extent[0]) and (x < grid.x_extent[1]) and 
#             (y > grid.y_extent[0]) and (y < grid.y_extent[1])):
#             print(x,y,lon,lat) #,footprint['tb'],old_footprint_list_15[i]['tb']
#             facecolor='blue'
#             ax = grid.plot_ellipse(x,y,
#                  major/250.0,minor/250.0,
#                  azim,
#                  color='black',alpha = 1.0,fig=fig0,facecolor=facecolor,ax=ax)

ax = grid.plot_ellipse(1.0,1.0,
                  major/250.0,major/250.0,
                  0.0,
                  color='black',alpha = 1.0,fig=fig0,facecolor='black',ax=ax)
ax.set_xlabel('EW Distance (km)')
ax.set_ylabel('NS Distance (km)')
ax.set_title('SSMIS, F18')
ax.text(1.4,0.9,'Desired Location')
ax.text(2.0,3.9,'Closest Location')

path = Path(f'L:/access/resampling/{sensor_name}/f{ksat}/location_plots')

png_file = path / 'footprint_locs_schematic.zoom.png'
fig0.savefig(png_file)
plt.show()
print