
import numpy as np
import matplotlib.pyplot as plt
import struct

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

def read_gain_array(sat_name='AMSR2',beamwidth=30,kscan=1,fov=1,gain_type='source',verbose = False):

    gain_array = np.zeros((1602,2402),dtype = 'float64')

    if gain_type == 'target':
        path = f"L:access/resampling/{sat_name}/target_gains/circular_{beamwidth:02d}km/"
        file_name = f"{path}s{kscan:02d}c{fov:04d}.dat"

    if gain_type == 'target_v2':
        path = f"L:access/resampling/{sat_name}/target_gains_v2/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'source':
        path = f"L:access/resampling/{sat_name}/source_gains/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'source_v2':
        path = f"L:access/resampling/{sat_name}/source_gains_v2/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'source_v3':
        path = f"L:access/resampling/{sat_name}/source_gains_v3/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'resample':
        path = f"L:access/resampling/{sat_name}/resample_gains/circular_{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:04d}.dat"

    elif gain_type == 'resample_v2':
        path = f"L:access/resampling/{sat_name}/resample_gains_v2/circular_{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:04d}.dat"

    else:
        raise ValueError(f'gain type {gain_type} makes no sense')


    print(file_name)
    num_lines_read = 0
    gain_tot = 0.0
    max_gain = 0.0

    with open(file_name,'rb') as f:
        for i in range(0,100000):
            try:
                x = f.read(12)

                ilat,ilon,gain = struct.unpack('<hhd',x)
              
                #if num_lines_read < 3:
                #    print(ilat,ilon,gain)
                num_lines_read += 1
                gain_array[ilat,ilon]  =gain
                gain_tot += gain
                if gain > max_gain:
                    max_gain = gain
                    ilat_max = ilat
                    ilon_max = ilon
                
            except:
                continue
    if verbose:
        print(f'Num Lines read = {num_lines_read}')
        print(f'Total Gain: {gain_tot:8.3f}')
        print()
        print(ilat_max,ilon_max)

    return gain_array,ilat_max,ilon_max


if __name__ == '__main__':
    plt.rcParams.update({'font.size': 16})
    target_lats,target_lons = read_location_file(sat_name='AMSR2',beamwidth=30)
    
    sat_name = 'AMSR2'
    beamwidth = 30
    band_num = 4
    
    plot_figs = True
    for kscan in [1]:
        #for fov in [1,5,11,20,21,243]:
        for fov in [1,2]:
            kscan_source =int((1+kscan)/2+14)
            kscan_source = 1
            fov_source = int((fov+1)/2)

            target_lat = target_lats[kscan-1,fov-1]
            target_lon = target_lons[kscan-1,fov-1]

            #These are wrong -- need to do projection....
            #target_ilon = 1 + np.floor((target_lon + 12.0)*100.0).astype('int')
            #target_ilat = 1 + np.floor((target_lat + 8.0)*100.0).astype('int')

            source_gain,ilat_max,ilon_max = read_gain_array(sat_name=sat_name,
                                                            beamwidth=beamwidth,
                                                            kscan=kscan_source,
                                                            fov = fov_source,
                                                            gain_type='source_v2',
                                                            verbose = False)
            source_gain3,ilat_max3,ilon_max3 = read_gain_array(sat_name=sat_name,
                                                            beamwidth=beamwidth,
                                                            kscan=kscan_source,
                                                            fov = fov_source,
                                                            gain_type='source_v3',
                                                            verbose = False)

            print(np.sum(source_gain))
            print(np.sum(source_gain3))                                            
            source_gain = source_gain/np.sum(source_gain)
            source_gain3 = source_gain3/np.sum(source_gain3)       


            l = 50
            source_sub =   source_gain[ilat_max-l:ilat_max+l+1,ilon_max-l:ilon_max+l+1]
            source_sub3 = source_gain3[ilat_max-l:ilat_max+l+1,ilon_max-l:ilon_max+l+1]
            
            max = np.nanmax(source_sub)
            source_sub = source_sub/max
            source_sub3 = source_sub3/max

            diff = source_sub3-source_sub

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
                target_fig,target_ax = plot_2d_array(source_sub3, xvals, yvals,zrange=(-1.0,1.0),
                        title=f'Source v3 Pattern, FOV = {fov:03d}', xtitle='EW distance (km)', ytitle='NS distance (km),',cmap='BrBG',plt_colorbar = True,fig=fig,subplot=222)
                target_ax.plot([0.0,0.0],[-l,l])
                target_ax.plot([-l,l],[0.0,0.0])
                # resample_fig,resample_ax = plot_2d_array(resample_sub, xvals, yvals,zrange=(-1.0,1.0), 
                #         title=f'Resampled Pattern, FOV = {fov:03d}', xtitle='EW distance (km)', ytitle='NS distance (km),',cmap='BrBG',plt_colorbar = True,fig=fig,subplot=223)
                # resample_ax.plot([0.0,0.0],[-l,l])
                # resample_ax.plot([-l,l],[0.0,0.0])
                diff_fig,diff_ax = plot_2d_array(diff, xvals, yvals,zrange=(-0.05,0.05), 
                        title=f'Source - Source V3 Pattern, FOV = {fov:03d}', xtitle='EW distance (km)', ytitle='NS distance (km),',cmap='BrBG',plt_colorbar = True,fig=fig,subplot=224)
                diff_ax.plot([0.0,0.0],[-l,l])
                diff_ax.plot([-l,l],[0.0,0.0])
                fig.tight_layout(h_pad=2)
                plt.show()
                # png_path = 'L:/access/resampling/AMSR2/resample_gains_v2/circular_30km/plots/'
                # png_file = f'{png_path}footprint_compare_band_{band_num:02d}_s{kscan:02d}_c{fov:03d}.png'
                # fig.savefig(png_file)
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



