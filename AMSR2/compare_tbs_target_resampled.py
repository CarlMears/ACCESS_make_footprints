


import sys
sys.path.append("C:/python_packages_cam/local_earth_grid")
sys.path.append("C:/job_access/python/land_water")
import numpy as np
import xarray as xr
import localearthgrid
from hansen_land_fraction import get_hansen_land_fraction_local
from read_gain_array import read_gain_array

import matplotlib.pyplot as plt
from numpy.random import RandomState


sat_name = 'AMSR2'
beamwidth = 30
band_num = 4

kscan = 1
fov = 243


target_gain,ilat_max_target,ilon_max_target = read_gain_array(sat_name=sat_name,
                                                beamwidth=beamwidth,
                                                band_num=band_num,
                                                kscan=kscan,
                                                fov=fov,
                                                gain_type='target_v2',
                                                verbose = True)
resample_gain,ilat_max2,ilon_max2 = read_gain_array(sat_name=sat_name,
                                                    beamwidth=beamwidth,
                                                    band_num=band_num,
                                                    kscan=kscan,
                                                    fov=fov,
                                                    gain_type='resample_v3',
                                                    verbose = True)  

target_gain = target_gain/np.sum(target_gain)
resample_gain = resample_gain/np.sum(resample_gain)
l = 50
lat_offset = 1
lon_offset = 0
target_sub    = target_gain[ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
resample_sub  = resample_gain[ilat_max_target-l+lat_offset:ilat_max_target+l+1+lat_offset,
                              ilon_max_target-l+lon_offset:ilon_max_target+l+1+lon_offset]

prng = RandomState(2730)
latlon_array = prng.rand(1000,2)
diff_list = []
for isamp in np.arange(0,1000):
    latlon = np.random.rand(2)
    lon = latlon_array[isamp,0]
    lat =latlon_array[isamp,1]
    lat = 30.0 + 40.0*lat
    lon = 230.0 + 40.0*lon
    print(f'lon = {lon}  lat = {lat}')
    lat_round_ten = 10.0*np.round(lat/10.0)
    lon_round_ten = 10.0*np.round(lon/10.0)
    if np.abs(lat-lat_round_ten) < 1.0:
        continue
    if np.abs(lon-lon_round_ten) < 1.0:
        continue

    grid = get_hansen_land_fraction_local(lon=lon,lat=lat,size_km = 50)
    grid.addlayer(name='Target')
    grid.setdata(target_sub,name='Target')
    grid.normalize(name='Target')

    grid.addlayer(name='Resamp')
    grid.setdata(resample_sub,name='Resamp')
    grid.normalize(name='Resamp')

    target_sum = np.sum(grid.multiply('Land Fraction','Target'))
    if target_sum < 0.1:
        continue
    resamp_sum = np.sum(grid.multiply('Land Fraction','Resamp'))

    plt.rcParams.update({'font.size': 18})
    from matplotlib import cm
    cmap = cm.winter
    diff = 150.0*abs(target_sum-resamp_sum)
    print(diff)
    if ((diff > 0.01) and (diff < 0.02)):
        ax = grid.contourplot(name='Land Fraction',vmin=0.0,vmax =1.0,
                         cbticks = [-1.0,-0.5,0.0,0.5,1.0],cmap = cmap,
                         units = 'Land Fraction',plt_contours=False,
                         xlabel='EW distance (km)',ylabel='NS distance (km)')

    # ax2 = grid.contourplot(name='Target',vmin=0.0,vmax = 1.0,
    #                     cbticks = [0.0,0.2,0.4,0.6,0.8,1.0],cmap = 'PuBu',
    #                     title = 'Target Pattern',plt_contours=False,
    #                     xlabel='EW distance (km)',ylabel='NS distance (km)')

    # ax3 = grid.contourplot(name='Resamp',vmin=0.0,vmax =1.0,
    #                     cbticks =  [0.0,0.2,0.4,0.6,0.8,1.0],cmap = 'PuBu',
    #                     title = 'Resampled Pattern',plt_contours=False,
    #                     xlabel='EW distance (km)',ylabel='NS distance (km)')


    print(f'Target: {target_sum:.5f} Resampled: {resamp_sum:.5f} Est Tb Diff: {150.0*(resamp_sum-target_sum):.5f}')
    diff_list.append(150.0*(resamp_sum-target_sum))
    plt.show()

diff_list = np.array(diff_list)
print(f'Mean: {np.mean(diff_list):.5f}   Std Dev: {np.std(diff_list):.5f}')
print(f'Max: {np.nanmax(diff_list):.5f}')
ax=plt.hist(diff_list,bins=100,range=[-1.0,1.0])

plt.show()
print()
