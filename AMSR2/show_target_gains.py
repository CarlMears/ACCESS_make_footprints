
import numpy as np
import matplotlib.pyplot as plt
import struct
import os
from read_gain_array import read_gain_array

def plot_2d_array(a, xvals, yvals,zrange=(0.0,1.0), title='', xtitle='', ytitle='',cmap='BrBG',plt_colorbar = True,zlabel=' ',fig=None,subplot=111):

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import copy

    X, Y = np.meshgrid(xvals, yvals)
    if fig is None:
        fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(subplot, title=title, xlabel=xtitle, ylabel=ytitle)
    im = ax.pcolormesh(X, Y, a, cmap=cmap, norm=colors.Normalize(vmin=zrange[0], vmax=zrange[1]))
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(20)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(16)

    if plt_colorbar:
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel(zlabel, fontsize=16)
    return fig,ax

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
    band_num = 4
    
    dlat = np.zeros((485))
    dlon = np.zeros((485))

    stddev_pattern_diff = np.zeros((485,2))
    tot_pattern_diff = np.zeros((485,2))

    plot_figs = True
    for kscan in [1]:
        #for fov in [1,5,11,20,21,243]:
        for fov in [1,5,11,20,21,242,243]:
            kscan_source =int((1+kscan)/2+14)
            fov_source = int((fov+1)/2)

            target_lat = target_lats[kscan-1,fov-1]
            target_lon = target_lons[kscan-1,fov-1]

            #These are wrong -- need to do projection....
            #target_ilon = 1 + np.floor((target_lon + 12.0)*100.0).astype('int')
            #target_ilat = 1 + np.floor((target_lat + 8.0)*100.0).astype('int')

            source_gain,ilat_max3,ilon_max3 = read_gain_array(sat_name=sat_name,
                                                            beamwidth=beamwidth,
                                                            band_num=band_num,
                                                            kscan=kscan_source,
                                                            fov = fov_source,
                                                            gain_type='source_v3',
                                                            verbose = False)
            target_gain,ilat_max,ilon_max = read_gain_array(sat_name=sat_name,
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
            resample_gain = resample_gain/np.sum(target_gain)

            ilat_max_source, ilon_max_source = np.unravel_index(source_gain.argmax(),source_gain.shape)
            ilat_max_target, ilon_max_target = np.unravel_index(target_gain.argmax(),target_gain.shape)
            ilat_max_resample, ilon_max_resample = np.unravel_index(resample_gain.argmax(),resample_gain.shape)

            lat_max_source = -8.0 + ilat_max_source*0.01
            lon_max_source = 0.0 - 12.0 + ilon_max_source*0.01

            lat_max_target = -8.0 + ilat_max_target*0.01
            lon_max_target = 0.0 - 12.0 + ilon_max_target*0.01

            lat_max_resample = -8.0 + ilat_max_resample*0.01
            lon_max_resample = 0.0 - 12.0 + ilon_max_resample*0.01


            print(f'Source Max Location = : {ilat_max_source:.3f}, {ilon_max_source:.3f}')       
            print(f'Source Max Location = : {lat_max_source:.3f}, {lon_max_source:.3f}')       
            print(f'Target Location = : {ilat_max_target:.3f}, {ilon_max_target:.3f}')
            print(f'Target Location = : {lat_max_target:.3f}, {lon_max_target:.3f}')
            print(f'Resample Location = : {ilat_max_resample:.3f}, {ilon_max_resample:.3f}')
            print(f'Resample Location = : {lat_max_resample:.3f}, {lon_max_resample:.3f}')


            l = 50
            source_sub = source_gain[ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
            target_sub = target_gain[ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
            resample_sub  = resample_gain [ilat_max_target-l:ilat_max_target+l+1,ilon_max_target-l:ilon_max_target+l+1]
            
            source_sub = source_sub/np.nanmax(source_sub)

            max = np.nanmax(target_sub)
            target_sub = target_sub/max
            resample_sub = resample_sub/max
            diff = resample_sub-target_sub

            if plot_figs:
                xvals = np.arange(-50,51)
                yvals = xvals

                # fig,ax = plt.subplots()
                # z = target_sub[50,:]
                # ax.plot(yvals,z)
                # z = target_sub[:,50]
                # ax.plot(xvals,z)
                # ax.plot([0.0,0.0],[0.0,np.nanmax(target_sub)])

                # fig,ax = plt.subplots()
                # z = source_sub[50,:]
                # ax.plot(yvals,z)
                # z = source_sub[:,50]
                # ax.plot(xvals,z)
                # ax.plot([0.0,0.0],[0.0,np.nanmax(source_sub)])
                fig = plt.figure(figsize=(14, 10))
                source_fig,source_ax = plot_2d_array(source_sub, xvals, yvals,zrange=(-1.0,1.0), 
                    title=f'Source Pattern, FOV = {fov_source:03d}', xtitle='EW distance (km)', ytitle='NS distance (km),',cmap='BrBG',plt_colorbar = True,fig=fig,subplot=221)
                source_ax.plot([0.0,0.0],[-l,l])
                source_ax.plot([-l,l],[0.0,0.0])
                target_fig,target_ax = plot_2d_array(target_sub, xvals, yvals,zrange=(-1.0,1.0),
                        title=f'Target Pattern, FOV = {fov:03d}', xtitle='EW distance (km)', ytitle='NS distance (km),',cmap='BrBG',plt_colorbar = True,fig=fig,subplot=222)
                target_ax.plot([0.0,0.0],[-l,l])
                target_ax.plot([-l,l],[0.0,0.0])
                resample_fig,resample_ax = plot_2d_array(resample_sub, xvals, yvals,zrange=(-1.0,1.0), 
                        title=f'Resampled Pattern, FOV = {fov:03d}', xtitle='EW distance (km)', ytitle='NS distance (km),',cmap='BrBG',plt_colorbar = True,fig=fig,subplot=223)
                resample_ax.plot([0.0,0.0],[-l,l])
                resample_ax.plot([-l,l],[0.0,0.0])
                diff_fig,diff_ax = plot_2d_array(diff, xvals, yvals,zrange=(-0.05,0.05), 
                        title=f'Resampled - Target Pattern, FOV = {fov:03d}', xtitle='EW distance (km)', ytitle='NS distance (km),',cmap='BrBG',plt_colorbar = True,fig=fig,subplot=224)
                diff_ax.plot([0.0,0.0],[-l,l])
                diff_ax.plot([-l,l],[0.0,0.0])
                fig.tight_layout(h_pad=2)

                png_path = 'L:/access/resampling/AMSR2/resample_gains_v2/circular_30km/plots/'
                os.makedirs(png_path,exist_ok = True)
                png_file = f'{png_path}footprint_compare_band_{band_num:02d}_s{kscan:02d}_c{fov:03d}.png'
                fig.savefig(png_file)
                plt.close(fig=fig)

            stddev_pattern_diff[fov-1,kscan-1] = np.std(diff)
            tot_pattern_diff[fov-1,kscan-1] = np.sum(diff)
            print(f'Sum of Pattern Difference = {np.sum(diff):.7f}')
            print(f'Std Dev of Pattern Difference = {np.std(diff):.7f}')
            print()

    fig,ax = plt.subplots()
    ax.plot(stddev_pattern_diff[:,0])
    ax.set_ylim([0.0,1.0e-6])

    fig2,ax2 = plt.subplots()
    ax2.plot(tot_pattern_diff[:,0])
    ax2.set_ylim([-0.02,0.02])
   
    print(np.nanmax(stddev_pattern_diff))

    plt.show()
    print()



