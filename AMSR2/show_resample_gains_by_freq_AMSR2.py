import sys

sys.path.append("M:/job_access/python/land_water/")
sys.path.append("M:/job_access/python/resample_wts/")
from rss_gridding.local_earth_grid import localearthgrid
from common.read_gain_array import read_gain_array
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from scipy.ndimage import shift

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

def read_debug_file(sat_name='AMSR2',ifreq=4):

    path = f"L:access/resampling/{sat_name}/source_gains/logs/"
    file_name = f"{path}find_source_gain_{sat_name}_band_{ifreq:02d}.debug"

    boresight_array = np.fromfile(file_name, dtype='float64',sep=' ')
    boresight_array = boresight_array.reshape((486,6))
    print()

    return boresight_array

def read_location_file(sat_name='AMSR2',beamwidth=30):

    path = f"L:access/resampling/{sat_name}/target_gains_v2/locs/"
    file_name = f"{path}find_target_gain_circular_{beamwidth:02d}km.loc"

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

    scan = np.reshape(scan,(2,485))
    fov = np.reshape(fov,(2,485))
    lats = np.reshape(lats,(2,485))
    lons = np.reshape(lons,(2,485))

    return lats,lons



if __name__ == '__main__':
    plt.rcParams.update({'font.size': 16})
    target_lats,target_lons = read_location_file(sat_name='AMSR2',beamwidth=30)
    
    sat_name = 'AMSR2'
    beamwidth = 30
    band_num = 1
    do_shift = True
    shift_vector = [1.41,1.41]

    panel_label_list = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n']

    
    dlat = np.zeros((485))
    dlon = np.zeros((485))

    stddev_pattern_diff = np.zeros((485,2,5))
    tot_pattern_diff = np.zeros((485,2,5))

    plot_figs = True
    freq_names = ['7GHz','11GHz','19GHz','24GHz','37GHz']
    #freq_names = ['11GHz','19GHz','24GHz','37GHz']
    
    for kscan in [1]:
        for fov in [1,5,11,20,21,31,41,151,243,244]:
        #for fov in np.arange(1,486):
            kscan_source =int((1+kscan)/2+14)
            fov_source = int((fov+1)/2)

            target_lat = target_lats[kscan-1,fov-1]
            target_lon = target_lons[kscan-1,fov-1]

            target_gain,ilat_max,ilon_max = read_gain_array(sat_name=sat_name,
                                                            beamwidth=beamwidth,
                                                            band_num=band_num,
                                                            kscan=kscan,
                                                            fov=fov,
                                                            gain_type='target_v2',
                                                            verbose = False)
            band_list = [1,2,3,4]
            resample_gain_array = np.zeros((1602,2402,len(band_list)))
            for iband,band_num in enumerate(band_list):
                resample_gain,ilat_max2,ilon_max2 = read_gain_array(sat_name=sat_name,
                                                                beamwidth=beamwidth,
                                                                band_num=band_num,
                                                                kscan=kscan,
                                                                fov=fov,
                                                                gain_type='resample_v3',
                                                                verbose = False) 
                resample_gain = resample_gain/np.sum(resample_gain)
                resample_gain_array[:,:,iband] = resample_gain
                 
            
            target_gain = target_gain/np.sum(target_gain)
            target_gain_shift = shift(target_gain,shift_vector)

            ilat_max_target, ilon_max_target = np.unravel_index(target_gain.argmax(),target_gain.shape)
            # ilat_max_resample, ilon_max_resample = np.unravel_index(resample_gain.argmax(),resample_gain.shape)

            # lat_max_source = -8.0 + ilat_max_source*0.01
            # lon_max_source = 0.0 - 12.0 + ilon_max_source*0.01

            # lat_max_target = -8.0 + ilat_max_target*0.01
            # lon_max_target = 0.0 - 12.0 + ilon_max_target*0.01

            # lat_max_resample = -8.0 + ilat_max_resample*0.01
            # lon_max_resample = 0.0 - 12.0 + ilon_max_resample*0.01


            # print(f'Source Max Location = : {ilat_max_source:.3f}, {ilon_max_source:.3f}')       
            # print(f'Source Max Location = : {lat_max_source:.3f}, {lon_max_source:.3f}')       
            # print(f'Target Location = : {ilat_max_target:.3f}, {ilon_max_target:.3f}')
            # print(f'Target Location = : {lat_max_target:.3f}, {lon_max_target:.3f}')
            # print(f'Resample Location = : {ilat_max_resample:.3f}, {ilon_max_resample:.3f}')
            # print(f'Resample Location = : {lat_max_resample:.3f}, {lon_max_resample:.3f}')


            l = 30
            #source_sub = source_gain[ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
            target_sub = target_gain[ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
            target_shift_sub = target_gain_shift[ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
            max = np.nanmax(target_sub)
            target_sub_to_plt = target_sub/max
            if do_shift:
                target_sub_to_plt = target_shift_sub/max
                

            if plot_figs:
                panel_index=0
                xvals = np.arange(-l,l+1)
                yvals = xvals

                fig = plt.figure(figsize=(25,10))

                target_fig,target_ax,target_im = plot_2d_array(target_sub_to_plt, xvals, yvals,zrange=(-1.0,1.0),
                                                     title=f'Target Pattern, FOV = {fov:03d}', 
                                                     xtitle='EW distance (km)', 
                                                     ytitle='NS distance (km)',
                                                     cmap='BrBG',
                                                     plt_colorbar = False,
                                                     fig=fig,
                                                     subplot=(2,5,1),
                                                     panel_label = panel_label_list[panel_index])
                panel_index += 1
                target_ax.plot([0.0,0.0],[-l,l])
                target_ax.plot([-l,l],[0.0,0.0])

            for iband,band_num in enumerate(band_list):
                plt_cbar = False
                if iband == len(band_list)-1:
                    plt_cbar=True
                resample_sub  = resample_gain_array[ilat_max_target-l:ilat_max_target+l+1,
                                                    ilon_max_target-l:ilon_max_target+l+1,
                                                    iband]
                resample_sub_to_plt = resample_sub/max
                if plot_figs:
                    resample_fig,resample_ax,resample_im = plot_2d_array(resample_sub_to_plt, xvals, yvals,zrange=(-1.0,1.0), 
                            title=f'FOV = {fov:03d}, {freq_names[band_num]}', 
                            xtitle='EW distance (km)', 
                            #ytitle='NS distance (km),',
                            cmap='BrBG',
                            plt_colorbar = False,
                            fig=fig,
                            subplot=(2,5,2+iband),
                            panel_label = panel_label_list[1+iband])
                    
                    resample_ax.plot([0.0,0.0],[-l,l])
                    resample_ax.plot([-l,l],[0.0,0.0])
                    resample_ax.yaxis.set_ticklabels([])

                if fov == 20:
                    print
                
                
                diff_for_stats = resample_sub/np.sum(resample_sub) - target_sub/np.sum(target_sub)
                diff_to_plot   = diff_for_stats/np.max(resample_sub)
                stddev_pattern_diff[fov-1,kscan-1,band_num] = np.std(diff_for_stats)
                if do_shift:
                    diff_for_stats = resample_sub/np.sum(resample_sub) - target_shift_sub/np.sum(target_shift_sub)
                    diff_to_plot   = diff_for_stats/np.max(resample_sub)
                    stddev_pattern_diff[fov-1,kscan-1,band_num] = np.std(diff_for_stats)

                print(f'Std Dev of Pattern Difference = {np.std(diff_for_stats):.7f}')
                print()

                if plot_figs:
                    diff_fig,diff_ax,diff_im = plot_2d_array(diff_to_plot, xvals, yvals,zrange=(-0.1,0.1), 
                                title=f'FOV = {fov:03d}, {freq_names[band_num]}', 
                                xtitle='EW distance (km)', 
                                #ytitle='NS distance (km),',
                                cmap='BrBG',
                                plt_colorbar = False,
                                fig=fig,
                                subplot=(2,5,7+iband),
                                panel_label = panel_label_list[5+iband])
                    panel_index += 1
                    diff_ax.plot([0.0,0.0],[-l,l])
                    diff_ax.plot([-l,l],[0.0,0.0])

            if plot_figs:
                cb_ax = fig.add_axes([0.92,0.57,0.015,0.3])
                cbar = fig.colorbar(target_im,cax=cb_ax,orientation='vertical')
                cbar.set_label('Amplitude',fontsize=12)
                cb_ax.tick_params(axis='y', labelsize=11)   

                cb_ax = fig.add_axes([0.92,0.12,0.015,0.3])
                cbar = fig.colorbar(diff_im,cax=cb_ax,orientation='vertical')
                cbar.set_label('Amplitude Difference',fontsize=12)
                cb_ax.tick_params(axis='y', labelsize=11)   
                #fig.tight_layout(h_pad=2)

                png_path = f'L:/access/resampling/AMSR2/resample_gains_v3/circular_{beamwidth:02d}km/plots/'
                os.makedirs(png_path,exist_ok = True)
                if do_shift:
                    png_file = f'{png_path}footprint_compare_band_all_bands_s{kscan:02d}_c{fov:03d}.shifted.png'
                    fig.savefig(png_file,bbox_inches='tight')
                    tif_file = f'{png_path}footprint_compare_band_all_bands_s{kscan:02d}_c{fov:03d}.shifted.tif'
                    fig.savefig(tif_file,bbox_inches='tight')
                else:
                    png_file = f'{png_path}footprint_compare_band_all_bands_s{kscan:02d}_c{fov:03d}.png'
                    fig.savefig(png_file,bbox_inches='tight')
                    tif_file = f'{png_path}footprint_compare_band_all_bands_s{kscan:02d}_c{fov:03d}.tif'
                    fig.savefig(tif_file,bbox_inches='tight')
                plt.close(fig=fig)



    # fig,ax = plt.subplots(figsize=(18,9))
    # for band_num in band_list:
    #     ax.plot(stddev_pattern_diff[:,0,band_num],label = freq_names[band_num])

    # ax.set_ylim([0.0,0.015])
    # ax.set_xlabel('Resampled FOV')
    # ax.set_ylabel('Std. Dev. (Resampled - Target)')
    # ax.legend()

    # png_path = Path(f'L:/access/resampling/AMSR2/resample_gains_v3/circular_{beamwidth:02d}km/plots/')
    # png_file = png_path / 'Std_dev_summary_by_fov.png'
    # fig.savefig(png_file)
    # tif_file = png_path.with_suffix('.tif')
    # fig.savefig(tif_file)
    plt.show()
    print()



