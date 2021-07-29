
import numpy as np
import matplotlib.pyplot as plt
import struct





def plot_2d_array(a, xvals, yvals,zrange=(0.0,1.0), title='', xtitle='', ytitle='',cmap='BrBG',plt_colorbar = True,zlabel=' '):

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import copy



    X, Y = np.meshgrid(xvals, yvals)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, title=title, xlabel=xtitle, ylabel=ytitle)
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

    path = f"L:access/resampling/{sat_name}/target_gains/locs/"
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

    elif gain_type == 'source':
        path = f"L:access/resampling/{sat_name}/source_gains/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'resample':
        path = f"L:access/resampling/{sat_name}/resample_gains/circular_{beamwidth:02d}km/band_{band_num:02d}/"
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

    sat_name = 'AMSR2'
    beamwidth = 30
    band_num = 2
    
    target_lats,target_lons = read_location_file(sat_name=sat_name,beamwidth=beamwidth)
    

    
    dlat = np.zeros((485))
    dlon = np.zeros((485))

    stddev_pattern_diff = np.zeros((485,2))
    tot_pattern_diff = np.zeros((485,2))
    lon_error_arr = np.zeros((485,2))
    lat_error_arr = np.zeros((485,2))
    l = 50

    dlat_arr = np.zeros((1+2*l,1+2*l))
    dlon_arr = np.zeros((1+2*l,1+2*l))

    for i in range(-l,l+1):
        dlat_arr[i+l,:] = i*0.01
        dlon_arr[:,i+l] = i*0.01



    plot_figs = False
    for kscan in [1,2]:
        for fov in range(1,486):
            kscan_source =kscan+14
            fov_source = int((fov+1)/2)

            target_lat = target_lats[kscan-1,fov-1]
            target_lon = target_lons[kscan-1,fov-1]

            target_ilon = 1 + np.floor((target_lon - 184.5 + 12.0)*100.0).astype('int')
            target_ilat = 1 + np.floor((target_lat + 8.0)*100.0).astype('int')

            target_gain,ilat_max,ilon_max = read_gain_array(sat_name=sat_name,beamwidth=beamwidth,kscan=kscan,fov=fov,gain_type='target',verbose = False)
            resample_gain,ilat_max2,ilon_max2 = read_gain_array(sat_name=sat_name,beamwidth=beamwidth,kscan=kscan,fov=fov,gain_type='resample',verbose = False)  
            
            target_gain = target_gain/np.sum(target_gain)
            resample_gain = resample_gain/np.sum(target_gain)

            ilat_max_target, ilon_max_target = np.unravel_index(target_gain.argmax(),target_gain.shape)
            ilat_max_resample, ilon_max_resample = np.unravel_index(resample_gain.argmax(),resample_gain.shape)

            lat_max_target = -8.0 + ilat_max_target*0.01
            lon_max_target = 184.5 - 12.0 + ilon_max_target*0.01

            lat_max_resample = -8.0 + ilat_max_resample*0.01
            lon_max_resample = 184.5 - 12.0 + ilon_max_resample*0.01


            print(f'Target Location = : {ilat_max_target:.3f}, {ilon_max_target:.3f}')
            print(f'Target Location = : {lat_max_target:.3f}, {lon_max_target:.3f}')
            print(f'Resample Location = : {ilat_max_resample:.3f}, {ilon_max_resample:.3f}')
            print(f'Resample Location = : {lat_max_resample:.3f}, {lon_max_resample:.3f}')


            l = 50
            target_sub = target_gain[target_ilat-l:target_ilat+l+1,target_ilon-l:target_ilon+l+1]
            resample_sub  = resample_gain [target_ilat-l:target_ilat+l+1,target_ilon-l:target_ilon+l+1]
            diff = resample_sub-target_sub

            if plot_figs:
                xvals = np.arange(-50,51)*0.01
                yvals = xvals

                fig,ax = plt.subplots()
                z = target_sub[50,:]
                ax.plot(yvals,z)
                z = target_sub[:,50]
                ax.plot(xvals,z)
                ax.plot([0.0,0.0],[0.0,np.nanmax(target_sub)])

                fig,ax = plt.subplots()
                z = source_sub[50,:]
                ax.plot(yvals,z)
                z = source_sub[:,50]
                ax.plot(xvals,z)
                ax.plot([0.0,0.0],[0.0,np.nanmax(source_sub)])




                source_fig,source_ax = plot_2d_array(source_sub, xvals, yvals,zrange=(-np.max(source_sub),np.max(source_sub)), 
                    title=f'Source Pattern, FOV = {fov_source:03d}', xtitle='Delta Lon.', ytitle='Delta Lat,',cmap='BrBG',plt_colorbar = True)
                source_ax.plot([0.0,0.0],[-0.5,0.5])
                source_ax.plot([-0.5,0.5],[0.0,0.0])
                target_fig,target_ax = plot_2d_array(target_sub, xvals, yvals,zrange=(-np.max(target_sub),np.max(target_sub)),
                        title=f'Target Pattern, FOV = {fov:03d}', xtitle='Delta Lon.', ytitle='Delta Lat,',cmap='BrBG',plt_colorbar = True)
                target_ax.plot([0.0,0.0],[-0.5,0.5])
                target_ax.plot([-0.5,0.5],[0.0,0.0])
                resample_fig,resample_ax = plot_2d_array(resample_sub, xvals, yvals,zrange=(-np.max(resample_sub),np.max(resample_sub)), 
                        title=f'Resampled Pattern, FOV = {fov:03d}', xtitle='Delta Lon.', ytitle='Delta Lat,',cmap='BrBG',plt_colorbar = True)
                resample_ax.plot([0.0,0.0],[-0.5,0.5])
                resample_ax.plot([-0.5,0.5],[0.0,0.0])
                diff_fig,diff_ax = plot_2d_array(diff, xvals, yvals,zrange=(-np.max(resample_sub)*0.1,np.max(resample_sub)*0.1), 
                        title=f'Resampled Pattern, FOV = {fov:03d}', xtitle='Delta Lon.', ytitle='Delta Lat,',cmap='BrBG',plt_colorbar = True)
                diff_ax.plot([0.0,0.0],[-0.5,0.5])
                diff_ax.plot([-0.5,0.5],[0.0,0.0])

            stddev_pattern_diff[fov-1,kscan-1] = np.std(diff)
            tot_pattern_diff[fov-1,kscan-1] = np.sum(diff)
            print(f'Sum of Pattern Difference = {np.sum(diff):.7f}')
            print(f'Std Dev of Pattern Difference = {np.std(diff):.7f}')

            lat_error_arr[fov-1,kscan-1] = np.sum(dlat_arr*resample_sub)
            lon_error_arr[fov-1,kscan-1] = np.sum(dlon_arr*resample_sub)
            print(f'Latitude Error = {lat_error_arr[fov-1,kscan-1]:.7f}')
            print(f'Longitude Error = {lon_error_arr[fov-1,kscan-1]:.7f}')

            print()

    plt_path = f'L:/access/resampling/{sat_name}/resample_weights/circular_{beamwidth:02d}km/plots/'

    fig,ax = plt.subplots()
    ax.plot(stddev_pattern_diff[:,0])
    ax.plot(stddev_pattern_diff[:,1])
    ax.set_ylim([0.0,3.0*np.mean(stddev_pattern_diff)])
    ax.set_xlim([0,486])
    ax.set_xlabel('FOV')
    ax.set_ylabel('Std. Dev. of Pattern Difference')
    png_file = f'{plt_path}stddev_pattern_diff_{band_num:02d}.png'
    fig.savefig(png_file)

    fig2,ax2 = plt.subplots()
    ax2.plot(tot_pattern_diff[:,0])
    ax2.plot(tot_pattern_diff[:,1])
    ax2.set_ylim([-3.0*np.mean(tot_pattern_diff),3.0*np.mean(tot_pattern_diff)])
    ax2.set_xlim([0,486])
    ax2.set_xlabel('FOV')
    ax2.set_ylabel('Total Pattern Difference')
    png_file = f'{plt_path}total_diff_{band_num:02d}.png'
    fig2.savefig(png_file)

    fig3,ax3 = plt.subplots()
    ax3.plot(lat_error_arr[:,0])
    ax3.plot(lat_error_arr[:,1])
    ax3.set_ylim([-0.02,0.02])
    ax3.set_xlim([0,486])
    ax3.set_xlabel('FOV')
    ax3.set_ylabel('Latitude Position Error')
    png_file = f'{plt_path}lat_error_{band_num:02d}.png'
    fig3.savefig(png_file)

    fig4,ax4 = plt.subplots()
    ax4.plot(lon_error_arr[:,0])
    ax4.plot(lon_error_arr[:,1])
    ax4.set_ylim([-0.02,0.02])
    ax4.set_xlim([0,486])
    ax4.set_xlabel('FOV')
    ax4.set_ylabel('Longitude Position Error')
    png_file = f'{plt_path}lon_error_{band_num:02d}.png'
    fig4.savefig(png_file)

    plt.show()
    print()



