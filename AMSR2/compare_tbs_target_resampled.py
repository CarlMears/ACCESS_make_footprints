


import sys

sys.path.append("M:/job_access/python/land_water/")
import numpy as np
import xarray as xr
from rss_gridding.local_earth_grid import localearthgrid
from hansen_land_fraction import get_hansen_land_fraction_local
from read_gain_array import read_gain_array

import matplotlib.pyplot as plt
from numpy.random import RandomState
from numpy.random import default_rng

def circle(*,diameter):
    phi = np.arange(0,2.0*np.pi,0.001)
    x = diameter/2.0 * np.cos(phi)
    y = diameter/2.0 * np.sin(phi)
    return (x,y)

rng = default_rng()
scene_list = ['lakes','midwest','coastline','edges','gradient']
band_list = [1,2,3,4]
band_names = ['7GHz','11GHz','19GHz','24GHz','37GHz']

rms_diff_array = np.zeros((len(scene_list),len(band_list)),dtype = np.float32)
mean_diff_array = np.zeros((len(scene_list),len(band_list)),dtype = np.float32)

scene_index = 0
scene = scene_list[scene_index]
sat_name = 'AMSR2'
land_water_contrast = 100.0

kscan = 1
fov = 243

beamwidth = 30
l = beamwidth

normalize=False

for band_index,band_num in enumerate(band_list):
    for scene_index,scene in enumerate(scene_list):

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


        target_sub    = target_gain[ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
        resample_sub  = resample_gain[ilat_max_target-l:ilat_max_target+l+1,
                                    ilon_max_target-l:ilon_max_target+l+1]
        
        if normalize:
            target_sub = target_sub/np.sum(target_sub)
            resample_sub = resample_sub/np.sum(resample_sub)

        dlatlon_array = rng.random((10000,10))
        diff_list = []
        target_list = []

        size_km = l
        isamp = 0
        num_successful = 0
        num_to_do = 250

        isamp=0
        while num_successful < num_to_do:
            if scene in ['lakes','midwest','coastline']:
                if scene == 'lakes':
                    center_lat = 52.0
                    center_lon = 296.0
                elif scene == 'midwest':
                    center_lat = 45.2
                    center_lon = 262.0
                elif scene == 'coastline':
                    center_lat = 38.5
                    center_lon = 236.8
                else:
                    raise ValueError(f'Scene: {scene} not implemented')
                dlon = -1.0 + 2.0*dlatlon_array[isamp,0]
                dlat = -1.0 + 2.0*dlatlon_array[isamp,1]
                if scene == 'coastline':
                    dlat = dlat/2.0
                    dlon = dlon/2.0

                isamp += 1
                lat = center_lat  + dlat
                lon = center_lon + dlon

                lat_round_ten = 10.0*np.round(lat/10.0)
                lon_round_ten = 10.0*np.round(lon/10.0)
                if np.abs(lat-lat_round_ten) < 1.0:
                    continue
                if np.abs(lon-lon_round_ten) < 1.0:
                    continue
        
                grid = get_hansen_land_fraction_local(lon=lon,lat=lat,size_km = size_km,verbose=False)
            elif scene == 'edges':
                # make grid with edge
                lon = 180.0
                lat = 0.0
                grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                                delta_x = 1.,delta_y = 1.,
                                num_x = 1+2*size_km,num_y = 1+2*size_km,
                                name='Land Fraction')
                phi = 2.0*np.pi*dlatlon_array[isamp,0]
                
                if beamwidth == 70.0:
                    offsetx = -20.0 + 40.0*dlatlon_array[isamp,1]
                    offsety = -20.0 + 40.0*dlatlon_array[isamp,2]
                else:
                    offsetx = -12.0 + 24.0*dlatlon_array[isamp,1]
                    offsety = -12.0 + 24.0*dlatlon_array[isamp,2]
                isamp += 1
                dx = np.cos(phi)
                dy = np.sin(phi)
                edge = np.zeros((1+size_km*2,1+size_km*2))
                for ix in range(-size_km,size_km+1):
                    for iy in range(-size_km,size_km+1):
                        if ((dx*(ix-offsetx))+(dy*(iy-offsety))) > 0.0:
                            edge[ix+size_km,iy+size_km] = 1.0
                grid.setdata(edge,'Land Fraction')
            elif scene == 'gradient':
                # make grid with randomly oriented gradient
                lon = 180.0
                lat = 0.0
                grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                                delta_x = 1.,delta_y = 1.,
                                num_x = 1+2*size_km,num_y = 1+2*size_km,
                                name='Land Fraction')
                phi = 2.0*np.pi*dlatlon_array[isamp,0]
                dx = np.cos(phi)
                dy = np.sin(phi)
                isamp += 1
                edge = np.zeros((1+size_km*2,1+size_km*2))
                for ix in range(-size_km,size_km+1):
                    for iy in range(-size_km,size_km+1):
                        dot = dx*ix +dy*iy
                        edge[ix+size_km,iy+size_km] = 0.5 + 0.01*dot
                grid.setdata(edge,'Land Fraction')
            else:
                raise ValueError(f'Scene: {scene} not implemented')
            
            grid.addlayer(name='Target')
            grid.setdata(target_sub,name='Target')
            grid.normalize(name='Target')

            grid.addlayer(name='Resamp')
            grid.setdata(resample_sub,name='Resamp')
            grid.normalize(name='Resamp')

            target_sum = np.sum(grid.multiply('Land Fraction','Target'))
            if scene in ['coastline','edges']:
                if target_sum < 0.15:
                    continue
                if target_sum > 0.85:
                    continue

            resamp_sum = np.sum(grid.multiply('Land Fraction','Resamp'))

            out_str = f'{scene}:{band_names[band_num]}: num = {num_successful} isamp= {isamp} lon = {lon:.2f}  lat = {lat:.2f}'
            out_str = f'{out_str} Target: {target_sum:.4f}' 
            out_str = f'{out_str} Resampled: {resamp_sum:.4f}'
            out_str = f'{out_str} Est Tb Diff: {land_water_contrast*(resamp_sum-target_sum):.5f}'
            print(out_str)

            diff_list.append(land_water_contrast*(resamp_sum-target_sum))
            target_list.append(target_sum)
            num_successful += 1
            
        target_list = np.array(target_list)
        diff_list = np.array(diff_list)

        print(f'Mean LF = {np.mean(target_list):.5f} +/- {np.std(target_list):.5f}')
        print(f'Range LF = {np.min(target_list):.5f} to {np.max(target_list):.5f}')

        print(f'Mean: {np.mean(diff_list):.5f}   RMS: {np.sqrt(np.mean(np.square(diff_list))):.5f}')
        print(f'Max: {np.nanmax(diff_list):.5f}')
        #ax=plt.hist(diff_list,bins=100,range=[-2.0,2.0])

        rms_diff_array[scene_index,band_index] = np.sqrt(np.mean(np.square(diff_list)))
        mean_diff_array[scene_index,band_index] = np.mean(diff_list)

    plt.show()
    print('RMS Difference')
    out_str = f'{band_names[band_num]}'
    for scene_index,scene in enumerate(scene_list):
        out_str = f'{out_str} {rms_diff_array[scene_index,band_index]:.5f}'
    print(out_str)
    print
    print('Mean Difference')
    out_str = f'{band_names[band_num]}'
    for scene_index,scene in enumerate(scene_list):
        out_str = f'{out_str} {mean_diff_array[scene_index,band_index]:.5f}'
    print(out_str)
    print
print()
print('---------------------------------------------------')
print(f'sat_name = {sat_name}')
print(f'land_water_contrast = {land_water_contrast}')
print(f'kscan = {kscan}')
print(f'fov = {fov}')
print(f'beamwidth = {beamwidth}')
print()
print('RMS Difference')
print(scene_list)
for band_index,band_num in enumerate(band_list):
    out_str = f'{band_names[band_num]}'
    for scene_index,scene in enumerate(scene_list):
        out_str = f'{out_str}   {rms_diff_array[scene_index,band_index]:.5f}'
    print(out_str)

print()
print('Mean Difference')
print(scene_list)
for band_index,band_num in enumerate(band_list):
    out_str = f'{band_names[band_num]}'
    for scene_index,scene in enumerate(scene_list):
        out_str = f'{out_str}   {mean_diff_array[scene_index,band_index]:.5f}'
    print(out_str)