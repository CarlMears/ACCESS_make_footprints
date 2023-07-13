
import sys
from pathlib import Path

sys.path.append("M:/job_access/python/land_water/")
sys.path.append("M:/job_access/python/resample_wts/")
import numpy as np
import xarray as xr
from rss_gridding.local_earth_grid import localearthgrid
from hansen_land_fraction import get_hansen_land_fraction_local
from common.read_gain_array import read_gain_array

import matplotlib.pyplot as plt
from numpy.random import RandomState
from numpy.random import default_rng

from common.random_land_mask import SomewhatRandomLandWater
from common.make_random_scene import MakeRandomScene

def circle(*,diameter):
    phi = np.arange(0,2.0*np.pi,0.001)
    x = diameter/2.0 * np.cos(phi)
    y = diameter/2.0 * np.sin(phi)
    return (x,y)


# scene_gen =SomewhatRandomLandWater()
scene_list = ['lakes','midwest','coastline','edges','gradient']
band_list = [1,2,4]
band_names = ['7GHz','11GHz','19GHz','24GHz','37GHz']
num_fovs=253

rms_diff_array = np.zeros((len(scene_list),len(band_list),num_fovs),dtype = np.float32)
mean_diff_array = np.zeros((len(scene_list),len(band_list),num_fovs),dtype = np.float32)

scene_index = 0
scene = scene_list[scene_index]
sat_name = 'AMSR2'
land_water_contrast = 100.0

kscan = 1
fov = 243


beamwidth = 30
subgrid_size = 50

num_to_do = 100

normalize=False
for scene_index,scene in enumerate(scene_list):
    scene_gen = MakeRandomScene(scene=scene)
    grids = []
    print(f'Loading {num_to_do} random land-water grids')
    for i in range(num_to_do):
        grid = scene_gen.make_random_scene(isamp = i)
        grids.append(grid)

    fig,ax = plt.subplots()
    for band_index,band_num in enumerate(band_list):
        for fov in range(1,num_fovs+1):
            target_gain,ilat_max_target,ilon_max_target = read_gain_array(sat_name=sat_name,
                                                            beamwidth=beamwidth,
                                                            band_num=band_num,
                                                            kscan=kscan,
                                                            fov=fov,
                                                            gain_type='target_v2',
                                                            verbose = False)
            resample_gain,ilat_max2,ilon_max2 = read_gain_array(sat_name=sat_name,
                                                                beamwidth=beamwidth,
                                                                band_num=band_num,
                                                                kscan=kscan,
                                                                fov=fov,
                                                                gain_type='resample_v3',
                                                                verbose = False)  
            
            target_gain = target_gain/np.sum(target_gain)
            resample_gain = resample_gain/np.sum(resample_gain)


            target_sub    = target_gain[ilat_max_target-subgrid_size:ilat_max_target+subgrid_size+1,ilon_max_target-subgrid_size:ilon_max_target+subgrid_size+1]
            resample_sub  = resample_gain[ilat_max_target-subgrid_size:ilat_max_target+subgrid_size+1,
                                        ilon_max_target-subgrid_size:ilon_max_target+subgrid_size+1]
            
            if normalize:
                target_sub = target_sub/np.sum(target_sub)
                resample_sub = resample_sub/np.sum(resample_sub)

            diff_list = []
            target_list = []

            for grid_index in range(num_to_do):
                grid,lat,lon = grids[grid_index]
                
                grid.addlayer(name='Target')
                grid.setdata(target_sub,name='Target')
                grid.normalize(name='Target')

                grid.addlayer(name='Resamp')
                grid.setdata(resample_sub,name='Resamp')
                grid.normalize(name='Resamp')

                target_sum = np.sum(grid.multiply('Land Fraction','Target'))
                resamp_sum = np.sum(grid.multiply('Land Fraction','Resamp'))

                out_str = f'{scene}:{band_names[band_num]}: num = {grid_index} lon = {grid.center_lon:.2f}  lat = {grid.center_lat:.2f}'
                out_str = f'{out_str} Target: {target_sum:.4f}' 
                out_str = f'{out_str} Resampled: {resamp_sum:.4f}'
                out_str = f'{out_str} Est Tb Diff: {land_water_contrast*(resamp_sum-target_sum):.5f}'
                #print(out_str)

                diff_list.append(land_water_contrast*(resamp_sum-target_sum))
                target_list.append(target_sum)
            
                
            target_list = np.array(target_list)
            diff_list = np.array(diff_list)

            #print(f'Mean LF = {np.mean(target_list):.5f} +/- {np.std(target_list):.5f}')
            #print(f'Range LF = {np.min(target_list):.5f} to {np.max(target_list):.5f}')

            print(f'Band: {band_names[band_num]} Scene: {scene} FOV: {fov}, Mean: {np.mean(diff_list):.5f}   RMS: {np.sqrt(np.mean(np.square(diff_list))):.5f}')
            #print(f'Max: {np.nanmax(diff_list):.5f}')
            #ax=plt.hist(diff_list,bins=100,range=[-2.0,2.0])

            rms_diff_array[scene_index,band_index,fov-1] = np.sqrt(np.mean(np.square(diff_list)))
            mean_diff_array[scene_index,band_index,fov-1] = np.mean(diff_list)

        ax.plot(np.arange(1,486),rms_diff_array[scene_index,band_index,:],label=band_names[band_num])

    ax.set_ylim(0.0,1.0)
    ax.set_xlabel('FOV index')
    ax.set_ylabel('Footprint Shape Induced Error (K)')
    ax.set_title(scene)
    ax.legend()
    png_file = Path(f'M:/job_access/docs/ResamplingPaper/figures/Shape_Error_{scene}_{beamwidth}km.png')
    fig.savefig(png_file)
    tif_file = png_file.with_suffix('.tif')
    fig.savefig(tif_file)

rms_error = xr.Dataset(
                data_vars=dict(
                    rms_error=(['scenes','bands','fovs'], rms_diff_array),
                    mean_error=(['scenes','bands','fovs'], mean_diff_array),
                ),
                coords=dict(
                    scenes=(["scenes"],np.arange(0,len(scene_list))),
                    bands=(["bands"],np.array(band_list)),
                    fov=(["fovs"],np.arange(0,num_fovs)),
                ),
                attrs=dict(
                    scene_names=scene_list,
                    beamwidth=beamwidth
                )

            )
nc_file = Path(f'M:/job_access/docs/ResamplingPaper/figures/Shape_Error_{beamwidth}km.nc')
rms_error.to_netcdf(nc_file)
plt.show()
print()