
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt

sensor_name = 'AMSR2'
ksat = 34
beamwidth = 30
kscan = 1
fov = 243

plot_path = Path('M:/job_access/python/resample_wts/plots/interpolation_study')
for band in [0,1,2,3,4]:
    fig1,ax1 = plt.subplots()
    for ib,beamwidth in enumerate([30,70]):
        if beamwidth < 60:
            band_offset = -1
            if band == 0:
                continue
        else:
            band_offset = 0
        nc_file = plot_path / f'rms_diff_table_{sensor_name}_{ksat}_{beamwidth}km_fov_{fov}_kscan_{kscan}_xr.nc'
        rms_diffs_xr = xr.open_dataset(nc_file)

        rms_diffs_one_band = rms_diffs_xr['rms_diff'].values[band+band_offset,:,:]
        scene_list = list(rms_diffs_xr['scene_names'].values)
        band_list = list(rms_diffs_xr['band_names'].values)

        x_offset_bw = (ib-1)*0.2
        bar_list = []
        bar = ax1.bar(np.arange(5)+x_offset_bw,rms_diffs_one_band[:,0],width=0.05,color='teal')
        bar_list.append(bar)
        bar = ax1.bar(np.arange(5)+x_offset_bw+0.05,rms_diffs_one_band[:,1],width=0.05,color='cyan')
        bar_list.append(bar)
        bar = ax1.bar(np.arange(5)+x_offset_bw+0.1,rms_diffs_one_band[:,2],width=0.05,color='salmon')
        bar_list.append(bar)
        for iscene in range(len(scene_list)):
            ax1.text(iscene+x_offset_bw,rms_diffs_one_band[iscene,2]*1.3,f'{beamwidth}km',rotation=90)
        
    scene_list.insert(0,'')
    scene_list = ['','lakes','farmland','coastline','edges','gradient']
    ax1.set_yscale('log')
    ax1.set_ylim((0.01,10))
    ax1.set_ylabel('Brightness Temperature Error (K)')
    ax1.set_xticklabels(scene_list)
    ax1.text(0.0,7.5,f'{band_list[band]}')
    ax1.legend(bar_list,['Shape Only','Interpolated','Closest'],loc='lower left',framealpha=1.0)

    png_file = plot_path / f'Resampling_Error_Summary_{sensor_name}_F{ksat}_{band_list[band]}.png'
    fig1.savefig(png_file)

plt.show()
print




        