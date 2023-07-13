import sys

sys.path.append("M:/job_access/python/land_water/")
sys.path.append("M:/job_access/python/resample_wts/")
from rss_gridding.local_earth_grid import LocalEarthGrid
from common.read_gain_array import read_gain_array
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from scipy.ndimage import shift

from SSMI.read_ssmi_location_file import read_location_file

def plot_2d_array(a, 
                  xvals, 
                  yvals,
                  zrange=(0.0,1.0), 
                  title='', 
                  xtitle='', 
                  ytitle='',
                  cmap='BrBG',
                  plt_colorbar = True,
                  zlabel=' ',
                  fig=None,
                  subplot=111,
                  panel_label=None,
                  panel_label_loc=[0.07,0.9]):

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import copy

    X, Y = np.meshgrid(xvals, yvals)
    if fig is None:
        fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(subplot[0],subplot[1],subplot[2], title=title, xlabel=xtitle, ylabel=ytitle)
    ax.set_aspect('equal')
    im = ax.pcolormesh(X, Y, a, cmap=cmap, norm=colors.Normalize(vmin=zrange[0], vmax=zrange[1]))
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(14)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(14)

    if plt_colorbar:
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel(zlabel, fontsize=12)

    if panel_label is not None:
        plt.text(panel_label_loc[0],
                 panel_label_loc[1],
                 panel_label,
                 transform=ax.transAxes,
                 fontsize=16)

    return fig,ax,im


if __name__ == '__main__':
    plt.rcParams.update({'font.size': 16})
    sat_name = 'SSMI'
    ksat=13
    beamwidth = 70
    freq_num = 0
    target_lats,target_lons = read_location_file(sat_name=sat_name,
                                                 ksat=ksat,
                                                 beamwidth=beamwidth,
                                                 freq=freq_num,
                                                 footprint_type='target')

    grid = LocalEarthGrid(center_lon = 0.0,center_lat = 0.0,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 2401,num_y = 1601,
                          name='Single')
    #panel_label_list = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n']
    
    num_scans = target_lats.shape[0]
    num_fovs = target_lats.shape[1]
    ix_diff = np.zeros((num_scans,num_fovs))
    iy_diff = np.zeros((num_scans,num_fovs))
    for scan in np.arange(0,num_scans):
        for fov in np.arange(0,num_fovs):
            target_gain,iy_max_target,ix_max_target = read_gain_array(sat_name=sat_name,
                                                                ksat=ksat,
                                                                beamwidth=beamwidth,
                                                                band_num=freq_num,
                                                                kscan=scan+1,
                                                                fov=fov+1,
                                                                gain_type='target',
                                                                verbose = False)
            
            z = grid.projection(target_lons[scan,fov],target_lats[scan,fov])
            ix_proj = z[0]/1000.0 + 1200.0
            iy_proj = z[1]/1000.0 + 800.0

            print(f'{scan+1:02d},{fov+1:03d} iy: {iy_proj:.2f}, {iy_max_target:.2f}, ix:{ix_proj:.2f}, {ix_max_target:.2f}')
            ix_diff[scan,fov] = ix_proj - ix_max_target
            iy_diff[scan,fov] = iy_proj - iy_max_target

    print('ix_stats',np.nanmin(ix_diff),np.nanmax(ix_diff),np.nanmean(ix_diff),np.nanstd(ix_diff))
    print('iy_stats',np.nanmin(iy_diff),np.nanmax(iy_diff),np.nanmean(iy_diff),np.nanstd(iy_diff))
    