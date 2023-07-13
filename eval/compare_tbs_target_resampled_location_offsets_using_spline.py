


import sys

sys.path.append("M:/job_access/python/land_water/")
sys.path.append("M:/job_access/python/resample_wts/")

import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
# from hansen_land_fraction import get_hansen_land_fraction_local
from common.read_gain_array import read_gain_array
from rss_gridding.local_earth_grid import LocalEarthGrid
from common.random_land_mask import SomewhatRandomLandWater

import matplotlib.pyplot as plt
from matplotlib import ticker

from numpy.random import RandomState
from numpy.random import default_rng

from scipy.ndimage import shift

from rss_gridding.quadrilateral_interpolation import InterpBilinearQuadrilateral

# orbit = 2
# tdr_path = 'L:/access/resampling/AMSR2/tdr/'
# tdr_file = f'{tdr_path}r{orbit:05d}.tdr.extra_fov_1_scan_1.h5'

def plot_sub_array(array_to_plot,name='Target'):
    
        sz = array_to_plot.shape
        num_x = sz[0]
        num_y = sz[1]

        grid = LocalEarthGrid(center_lon = target_lons[kscan-1,fov-1],center_lat = target_lats[kscan-1,fov-1],
                                     delta_x = 1.,delta_y = 1.,
                                    num_x = num_x,num_y = num_y,
                                    name=name)
        grid.setdata(array_to_plot,name)
        grid.normalize(name=name)
        maxval = np.nanmax(array_to_plot)
        minval = -maxval
        fig,ax,im = grid.contourplot(name=name,title =name,vmin=-maxval,vmax=maxval,scale=False)

        return fig,ax,im

def rms(x):
    return np.sqrt(np.mean(x*x))

# def find_AMSR2_target_locations():

#     print(tdr_file)

#     d = xr.open_dataset(tdr_file)
#     lats = d.latitude.values
#     lons = d.longitude.values
#     azim = d.azimuth.values

#     #for this orbit, the scan places the average lat of all 29x243 scans
#     #very close to the equator

#     scan0 = 1857
#     lats = lats[scan0:scan0+2,:]
#     lons = lons[scan0:scan0+2,:]-lons[scan0,242]
#     azim = azim[scan0:scan0+2,:]

#     return lats,lons,azim

def load_equatorial_target_locs(sensor_name = 'AMSR2',ksat=18):

    if sensor_name == 'AMSR2':
        path_to_loc_file = Path(f'L:/access/resampling/AMSR2/target_gains_v2/locs/find_target_gain_circular_30km.loc')
        num_scans = 2
        num_fovs = 485
    elif sensor_name == 'SSMIS':
        path_to_loc_file = \
            Path(f'L:/access/resampling/SSMIS/f{ksat:02d}/target_gains/locs/find_target_gain_f{ksat:02d}_circular.loc')
        num_scans = 2
        num_fovs = 446
    else:
        raise ValueError(f'sensor_name = {sensor_name} invalid')

    lats = np.zeros((num_scans,num_fovs),dtype=np.float32)
    lons = np.zeros((num_scans,num_fovs),dtype=np.float32)
    with open(path_to_loc_file,'r') as f:
        lines = f.readlines()
        for line in lines:
            items = line.split()
            scan = int(items[0])
            fov = int(items[1])-1
            lat = float(items[2])
            lon = float(items[3])

            lats[scan,fov] = lat
            lons[scan,fov] = lon

    return lats,lons

def circle(*,diameter):
    phi = np.arange(0,2.0*np.pi,0.001)
    x = diameter/2.0 * np.cos(phi)
    y = diameter/2.0 * np.sin(phi)
    return (x,y)


if __name__ == '__main__':
    #get the target locations
    #target_lats,target_lons,azim = find_AMSR2_target_locations() 

    sensor_name = 'SSMI:'
    ksat = 13
    beamwidth = 70
    plot_example_footprints = True # set true to remake Figure 7
    plot_path = Path('L:/access/resampling/SSMI/f13/eval/interpolation_study')


    target_lats,target_lons = load_equatorial_target_locs(sensor_name = sensor_name,ksat=ksat)

    #set up the local grid -- used to find distances in km in swath-relative coordinates
    grid = LocalEarthGrid(center_lon = 0.0,center_lat = 0.0,
                            delta_x = 1.,delta_y = 1.,
                            num_x = 2401,num_y = 1601,
                            name='Single')

    rng = default_rng()
    scene_list = ['lakes','midwest','coastline','edges','gradient']
    
    if sensor_name == 'AMSR2':
        if beamwidth == 30:
            band_name = ['7GHz','11GHz','19GHz','24GHz','36GHz']
            band_list = [0,1,2,3]
            band_list = [0,1]  #Temporary for testing
            freq_index_offset = 1
        elif beamwidth == 70:
            band_name = ['7GHz','11GHz','19GHz','24GHz','36GHz']
            band_list = [0,1,2,3,4]
            freq_index_offset = 0
        else:
            raise ValueError(f'Beamwidth {beamwidth} not valid')

        kscan = 1
        center_fov = 243
    elif sensor_name == 'SSMIS':
        band_name = ['19GHz','22GHz','37GHz']
        band_list = [0,1,2]
        freq_index_offset = 2
        kscan = 1
        center_fov = 223
    elif sensor_name == 'SSMI':
        band_name = ['19GHz','22GHz','37GHz']
        band_list = [0,1,2]
        freq_index_offset = 2
        kscan = 1
        center_fov = 223
    else:
        raise ValueError(f'sensor name {sensor_name} invalid')

    land_water_contrast = 100.0
    gradient_magnitude = 0.5  #K/kilometer or 50K / 100 km
    num_to_do = 1000

    kscan = 1
    fov = 243
    xlocs = np.zeros((2,2))
    ylocs = np.zeros((2,2))

    for dscan in [0,1]:
        for dfov in [0,1]:
            xlocs[dscan,dfov],ylocs[dscan,dfov] = grid.projection(target_lons[kscan+dscan-1,fov-dfov-1],target_lats[kscan+dscan-1,fov-dfov-1])

    #convert to km
    xlocs = xlocs/1000.0
    ylocs = ylocs/1000.0

    x_locs_ordered = [xlocs[0,0],xlocs[0,1],xlocs[1,1],xlocs[1,0]]
    y_locs_ordered = [ylocs[0,0],ylocs[0,1],ylocs[1,1],ylocs[1,0]]

    #x_locs_ordered = [xlocs[0,0],xlocs[1,0],xlocs[1,1],xlocs[0,1]]
    #y_locs_ordered = [ylocs[0,0],ylocs[1,0],ylocs[1,1],ylocs[0,1]]
    z = np.zeros((4))
    qi = InterpBilinearQuadrilateral(x_locs_ordered,y_locs_ordered,z)
    
    center_x = np.sum(xlocs)/4.0
    center_y = np.sum(ylocs)/4.0

    # this is the size of the rectangle around the 
    y_extent = 2.0*np.abs(center_y-ylocs[0,0])
    x_extent = 2.0*np.abs(center_x-xlocs[0,0])

    size_km = beamwidth
    use_precomputed_locations = False
    subgrid_size = 140

    std_diff_table_ideal = np.zeros((len(band_list),len(scene_list)))
    std_diff_table_closest = np.zeros((len(band_list),len(scene_list)))
    std_diff_table_interp = np.zeros((len(band_list),len(scene_list)))

    rms_diff_table = np.zeros((len(band_list),len(scene_list),3))
    mean_diff_table = np.zeros((len(band_list),len(scene_list)))

    diff_list_ideal_all = np.full((len(band_list),len(scene_list),num_to_do),np.nan)
    diff_list_closest_all = np.full((len(band_list),len(scene_list),num_to_do),np.nan)
    diff_list_interp_all = np.full((len(band_list),len(scene_list),num_to_do),np.nan)
    target_list_all = np.full((len(band_list),len(scene_list),num_to_do),np.nan)
    distance_closest_all = np.full((len(band_list),len(scene_list),num_to_do),np.nan)

    for iband,band_num in enumerate(band_list):
        if iband == 0:
            continue
        if sensor_name == 'AMSR2':
            target_gain,ilat_max_target,ilon_max_target = read_gain_array(sat_name=sensor_name,
                                                    beamwidth=beamwidth,
                                                    band_num=band_num+freq_index_offset,
                                                    kscan=kscan,
                                                    fov=fov,
                                                    gain_type='target_v2',
                                                    verbose = True)
            #if beamwidth > 60:
            #    use_precomputed_locations=True
        elif sensor_name in ['SSMIS','SSMI']:
            target_gain,ilat_max_target,ilon_max_target = read_gain_array(sat_name=sensor_name,
                                                    ksat=ksat,
                                                    beamwidth=beamwidth,
                                                    band_num=band_num+freq_index_offset,
                                                    kscan=kscan,
                                                    fov=fov,
                                                    gain_type='target',
                                                    verbose = True)

        target_gain = target_gain/np.sum(target_gain)
        
        resample_sub_array = np.zeros((1+2*subgrid_size,1+2*subgrid_size,2,2))

        for dscan in [0,1]:
            for dfov in [0,1]:
                if sensor_name == 'AMSR2':
                    resample_gain,ilat_max2,ilon_max2 = read_gain_array(sat_name=sensor_name,
                                                        beamwidth=beamwidth,
                                                        band_num=band_num+freq_index_offset,
                                                        kscan=kscan+dscan,
                                                        fov=fov-dfov,
                                                        gain_type='resample_v3',
                                                        verbose = True)  
                elif sensor_name in ['SSMIS','SSMI']:
                    resample_gain,ilat_max2,ilon_max2 = read_gain_array(sat_name=sensor_name,
                                                        ksat=ksat,
                                                        beamwidth=beamwidth,
                                                        band_num=band_num+freq_index_offset,
                                                        kscan=kscan+dscan,
                                                        fov=fov-dfov,
                                                        gain_type='resample',
                                                        verbose = True) 
                else:
                    raise ValueError(f'sensor_name {sensor_name} invalid')

                resample_gain = resample_gain/np.sum(resample_gain)
                resample_sub  = resample_gain[ilat_max_target-subgrid_size:ilat_max_target+subgrid_size+1,
                                              ilon_max_target-subgrid_size:ilon_max_target+subgrid_size+1]
                title_str = f''
                resample_sub_array[:,:,dscan,dfov] = resample_sub

        target_sub = target_gain[ilat_max_target-subgrid_size:ilat_max_target+subgrid_size+1,
                                 ilon_max_target-subgrid_size:ilon_max_target+subgrid_size+1]

        for scene_index,scene in enumerate(scene_list):
            if scene in ['lakes','midwest','coastline']:
                scene_gen =SomewhatRandomLandWater(scene=scene)

            diff_list_ideal = []
            diff_list_closest = []
            diff_list_interp = []
            distance_list_closest = []

            print(f'Computing Scene = {scene}, Band = {band_name[band_num+freq_index_offset]}')
            
            num_good = 0
            num_random_samples=10000
            num_done = 0
            lats_all = []
            lons_all = []
            phi_all = []
            shift_x_all = []
            shift_y_all = []

            rng = np.random.default_rng()
            random_array = rng.random((num_random_samples,10))
            isamp = 0

            if use_precomputed_locations:
                    txt_file = plot_path / f'random_locations_{scene}_{60}km.txt'
                    df = pd.read_csv(txt_file,header=None)
                    lons = df.iloc[:,0].values
                    lats = df.iloc[:,1].values
                    xshft = df.iloc[:,2].values
                    yshft = df.iloc[:,3].values
                    if scene in ['edges','gradient']:
                        phis = df.iloc[:,4].values

            while num_done < num_to_do:  
                if scene in ['lakes','midwest','coastline']:
                    if use_precomputed_locations:
                        grid,lat,lon = scene_gen.get_somewhat_random_earth_scene(subgrid_size=subgrid_size,
                                                                                 lon_in=lons[isamp],
                                                                                 lat_in = lats[isamp])
                    else:
                        grid,lat,lon = scene_gen.get_somewhat_random_earth_scene(subgrid_size=subgrid_size)

                elif scene == 'edges':
                    # make grid with edge

                    # for the edge case, we want to make a grid with a single edge
                    # it needs to be at a random angle, and a random distance 
                    # from the center.  The distance is chosen to be within 1/2 of
                    # the beamwidth.  The angle is chosen randomly.

                    lon = 180.0
                    lat = 0.0
                    grid = LocalEarthGrid(center_lon = lon,center_lat = lat,
                                    delta_x = 1.,delta_y = 1.,
                                    num_x = 1+2*subgrid_size,num_y = 1+2*subgrid_size,
                                    name='Land Fraction')
                    
                    # chose a random angle for the edge
                    if use_precomputed_locations:
                            phi = phis[isamp]
                    else:
                        phi = 2.0*np.pi*random_array[isamp,0]

                    # chose a random offset from the center
                    offsetx = -0.5*beamwidth + beamwidth*random_array[isamp,1]
                    offsety = -0.5*beamwidth + beamwidth*random_array[isamp,2]
                    
                    
                    phi_all.append(phi)
                    dx = np.cos(phi)
                    dy = np.sin(phi)

                    edge = np.zeros((2*subgrid_size+1,2*subgrid_size+1))
                    for ix in range(-subgrid_size,subgrid_size+1):
                        for iy in range(-subgrid_size,subgrid_size+1):
                            if ((dx*(ix-offsetx))+(dy*(iy-offsety))) > 0.0:
                                edge[ix+subgrid_size,iy+subgrid_size] = 1.0
                    grid.setdata(edge,'Land Fraction')

                elif scene == 'gradient':
                    # make grid spatial gradient

                    # for the gradient case, the gradient is along a random angle
                    # translation is not required

                    lon = 180.0
                    lat = 0.0
                    grid = LocalEarthGrid(center_lon = lon,center_lat = lat,
                                    delta_x = 1.,delta_y = 1.,
                                    num_x = 1+2*subgrid_size,num_y = 1+2*subgrid_size,
                                    name='Land Fraction')

                    if use_precomputed_locations:
                        phi = phis[isamp]
                    else:
                        phi = 2.0*np.pi*random_array[isamp,0]
 
                    phi_all.append(phi)
                    dx = np.cos(phi)
                    dy = np.sin(phi)
                    grad = np.zeros((2*subgrid_size+1,2*subgrid_size+1))
                    for ix in range(-subgrid_size,subgrid_size+1):
                        for iy in range(-subgrid_size,subgrid_size+1):
                            dot = dx*ix +dy*iy
                            grad[ix+subgrid_size,iy+subgrid_size] = (gradient_magnitude/land_water_contrast)*dot
                    grid.setdata(grad,'Land Fraction')
                else:
                    raise ValueError(f'Scene: {scene} not implemented')

                dfov_vector = np.array([xlocs[0,1] -xlocs[0,0],ylocs[0,1] -ylocs[0,0]])
                dscan_vector = np.array([xlocs[1,0] -xlocs[0,0],ylocs[1,0] -ylocs[0,0]])
                
                if use_precomputed_locations:
                    shift_vector = [xshft[isamp],yshft[isamp],]
                else:
                    shift_vector = random_array[isamp,6]*dfov_vector + \
                                   random_array[isamp,7]*dscan_vector

                if plot_example_footprints:
                    shift_vector = [2.0,2.0]

                shifted_loc = [xlocs[0,0]+shift_vector[0],ylocs[0,0]+shift_vector[1]]

                isamp+=1
                shift_x_all.append(shift_vector[0])
                shift_y_all.append(shift_vector[1])

                # compute the shifted target pattern

                target_gain_shift = shift(target_gain,shift_vector[::-1])
                target_sub_shift = target_gain_shift[ilat_max_target-subgrid_size:ilat_max_target+subgrid_size+1,
                                                     ilon_max_target-subgrid_size:ilon_max_target+subgrid_size+1]
        
                # plot_sub_array(target_sub,name='Target')
                # plot_sub_array(target_sub_shift,name='Target')
                # print
                #find the interpolation weights for this shift
                try:
                    wts = qi.weights(xlocs[0,0]+shift_vector[0],ylocs[0,0]+shift_vector[1])
                except ValueError:
                    print(xlocs)
                    print(ylocs)
                    print(shift_vector)
                    print
                wts_array = np.zeros((2,2))
                wts_array[0,0] = wts[0]
                wts_array[0,1] = wts[1]
                wts_array[1,1] = wts[2]
                wts_array[1,0] = wts[3]

                closest_corner = np.unravel_index(np.argmax(wts_array),wts_array.shape)

                for dscan in [0,1]:
                    for dfov in [0,1]:
                        distance_to_shifted_loc = np.sqrt(np.square(xlocs[dscan,dfov]-shifted_loc[0]) +
                                                          np.square(ylocs[dscan,dfov]-shifted_loc[1]))
                        # print(f'{dscan},{dfov}: {distance_to_shifted_loc:.2f}, wt: {wts_array[dscan,dfov]:.3f}')
                
                # print(closest_corner)
                # print(f'Weight, closest corner {wts_array[closest_corner]:.3f}')
                distance_closest = np.sqrt(np.square(xlocs[closest_corner]-shifted_loc[0]) +
                       np.square(ylocs[closest_corner]-shifted_loc[1]))
                # print(f'Distance, closest corner {distance_closest:.3f}')
                        
                # print

                #find the interpolated resampled pattern
                resample_sub_interp = np.zeros_like(resample_sub)
                for dscan in [0,1]:
                    for dfov in [0,1]:
                        #resample_sub_interp = resample_sub_interp + resample_sub_array[:,:,dfov,dscan]*wts_array[dfov,dscan]
                        resample_sub_interp = resample_sub_interp + resample_sub_array[:,:,dscan,dfov]*wts_array[dscan,dfov]
                
                resample_sub_closest = resample_sub_array[:,:,closest_corner[0],closest_corner[1]]

                grid.addlayer(name='Target')
                grid.setdata(target_sub,name='Target')
                grid.normalize(name='Target')

                grid.addlayer(name='Target_Shift')
                grid.setdata(target_sub_shift,name='Target_Shift')
                grid.normalize(name='Target')

                grid.addlayer(name='Resamp')
                grid.setdata(resample_sub_array[:,:,0,0],name = 'Resamp')
                grid.normalize(name='Resamp')


                grid.addlayer(name='Resamp_Closest')
                grid.setdata(resample_sub_closest,name='Resamp_Closest')
                grid.normalize(name='Resamp_Closest')

                grid.addlayer(name='Resamp_Interp')
                grid.setdata(resample_sub_interp,name='Resamp_Interp')
                grid.normalize(name='Resamp_Interp')

                grid.addlayer(name='Diff_Ideal')
                grid.setdata(target_sub-resample_sub_array[:,:,0,0],name='Diff_Ideal')

                grid.addlayer(name='Diff_Shifted')
                grid.setdata(target_sub_shift-resample_sub_array[:,:,0,0],name='Diff_Shifted')

                grid.addlayer(name='Diff_Interp')
                grid.setdata(target_sub_shift-resample_sub_interp[:,:],name='Diff_Interp')

                if plot_example_footprints:

                    # This makes Figure 08 in the resampling paper

                    plt.rcParams.update({'font.size': 14})
                    plt.rcParams["xtick.direction"] = "in"
                    plt.rcParams["ytick.direction"] = "in"
                    fig,axs = plt.subplots(nrows=3,ncols=2,figsize=(9,12))
                    maxval = np.max(target_sub)*1000.  #going to scale the footprints later
                    cbticks = np.array([-0.001,-0.0005,0.0,0.0005,0.001])
                    fig,axs[0,0],_ = grid.contourplot(fig_in=fig,
                                                 ax_in = axs[0,0],
                                                 name='Target',
                                                 title ='Target',
                                                 panel_label = '(a)',
                                                 plt_contours=False,
                                                 plt_colorbar=False,
                                                 vmin=-maxval,
                                                 vmax=maxval,
                                                 scale=True,
                                                 scale_factor=0.001,
                                                 xlabel = '',
                                                 cbticks=cbticks)
                    
                    axs[0,0].plot([-50,50],[0,0],color='gray')
                    axs[0,0].plot([0,0],[-50,50],color='gray')

                    fig,ax,_ = grid.contourplot(fig_in=fig,
                                                 ax_in = axs[0,1],
                                                 name='Target_Shift',
                                                 title ='Target, Shifted',
                                                 panel_label = '(b)',
                                                 vmin=-maxval,
                                                 vmax=maxval,
                                                 scale=True,
                                                 scale_factor=0.001,
                                                 xlabel='',
                                                 ylabel='',
                                                 plt_contours=False,
                                                 plt_colorbar=False,
                                                 cbticks=cbticks)
                    axs[0,1].plot([-50,50],[0,0],color='gray')
                    axs[0,1].plot([0,0],[-50,50],color='gray')
                    axs[0,1].set_yticklabels(' ')

                    fig,ax,_ = grid.contourplot(fig_in=fig,
                                                 ax_in = axs[1,0],
                                                 name='Resamp',
                                                 title ='Resampled, Closest',
                                                 panel_label = '(c)',
                                                 vmin=-maxval,
                                                 vmax=maxval,
                                                 scale=True,
                                                 scale_factor=0.001,
                                                 plt_colorbar=False,
                                                 plt_contours=False,
                                                 xlabel='',
                                                 cbticks=cbticks)
                    axs[1,0].plot([-50,50],[0,0],color='gray')
                    axs[1,0].plot([0,0],[-50,50],color='gray')

                    fig,ax,im_foot = grid.contourplot(fig_in=fig,
                                                 ax_in = axs[1,1],
                                                 name='Resamp_Interp',
                                                 title ='Resampled, Interpolated',
                                                 panel_label = '(d)',
                                                 vmin=-maxval,
                                                 vmax=maxval,
                                                 scale=True,
                                                 scale_factor=0.001,
                                                 xlabel='',
                                                 ylabel='',
                                                 plt_colorbar=False,
                                                 plt_contours=False,
                                                 cbticks=cbticks)
                    axs[1,1].plot([-50,50],[0,0],color='gray')
                    axs[1,1].plot([0,0],[-50,50],color='gray')
                    axs[1,1].set_yticklabels(' ')
                    
                    # fig,ax,_ = grid.contourplot(fig_in=fig,
                    #                              ax_in = axs[2,0],
                    #                              name='Diff_Ideal',
                    #                              title ='Difference, Ideal',
                    #                              panel_label = '(e)',
                    #                              vmin=-maxval/10.,
                    #                              vmax=maxval/10.,
                    #                              scale=True,
                    #                              scale_factor=0.001,
                    #                              plt_colorbar=False,
                    #                              plt_contours=False,
                    #                              cbticks=cbticks/10.)
                    # axs[2,0].plot([-50,50],[0,0],color='gray')
                    # axs[2,0].plot([0,0],[-50,50],color='gray')
            
                    fig,ax,_ = grid.contourplot(fig_in=fig,
                                                 ax_in = axs[2,0],
                                                 name='Diff_Shifted',
                                                 title ='Difference, Closest',
                                                 panel_label = '(e)',
                                                 vmin=-maxval/10.,
                                                 vmax=maxval/10.,
                                                 scale=True,
                                                 scale_factor=0.001,
                                                 ylabel='',
                                                 plt_colorbar=False,
                                                 plt_contours=False,
                                                 cbticks=cbticks/10.)
                    axs[2,0].plot([-50,50],[0,0],color='gray')
                    axs[2,0].plot([0,0],[-50,50],color='gray')
                    
                    fig,ax,im_diff = grid.contourplot(fig_in=fig,
                                                 ax_in = axs[2,1],
                                                 name='Diff_Interp',
                                                 title ='Difference, Interpolated',
                                                 panel_label = '(f)',
                                                 vmin=-maxval/10.,
                                                 vmax=maxval/10.,
                                                 scale=True,
                                                 scale_factor=0.001,
                                                 ylabel='',
                                                 plt_colorbar=False,
                                                 plt_contours=False,
                                                 cbticks=cbticks/10.)
                    
                    axs[2,1].plot([-50,50],[0,0],color='gray')
                    axs[2,1].plot([0,0],[-50,50],color='gray')
                    axs[2,1].set_yticklabels(' ')

                    fig.subplots_adjust(bottom=0.1,top=0.95,
                                        left = 0.1,right=0.8,
                                        hspace=0.25,wspace=0.05)
                    cb_ax = fig.add_axes([0.85,0.5,0.025,0.4])
                    cbar = fig.colorbar(im_foot,cax=cb_ax,orientation='vertical')
                    cbar.set_label(f'Footprint Amplitude')
                    cb_ax.tick_params(axis='x', labelsize=11)
                    #cb_ax.ticklabel_format(scilimits=(0,0))

                    cb_ax2 = fig.add_axes([0.85,0.125,0.025,0.22])
                    cbar = fig.colorbar(im_diff,cax=cb_ax2,orientation='vertical')
                    cbar.set_label(f'Footprint Difference')
                    cb_ax2.tick_params(axis='x', labelsize=11)
                    
                    png_file = plot_path / f'footprints_combined_{beamwidth:03d}.png'
                    fig.savefig(png_file)

                    # plot_path_fig = Path('M:/job_access/docs/ResamplingPaper/figures/Figure_08_Location_Error_Explain/')
                    # png_file = plot_path_fig / 'Fig_08.png'
                    # fig.savefig(png_file)
                    # tif_file = png_file.with_suffix('.tif')
                    # fig.savefig(tif_file,dpi=300)


                    plt.rcParams.update({'font.size': 18})
                    plt.rcParams["xtick.direction"] = "out"
                    plt.rcParams["ytick.direction"] = "out"
                    plot_example_footprints = False
        
                target_sum = np.sum(grid.multiply('Land Fraction','Target'))
                target_sum_shift = np.sum(grid.multiply('Land Fraction','Target_Shift'))
                if not use_precomputed_locations:
                    # ensure that the target location is actually coastline
                    if scene in ['coastline','edges']:
                        if target_sum_shift < 0.15:
                            continue
                        if target_sum_shift > 0.85:
                            continue

                lats_all.append(lat)
                lons_all.append(lon)
                resamp_sum = np.sum(grid.multiply('Land Fraction','Resamp'))
                resamp_sum_closest = np.sum(grid.multiply('Land Fraction','Resamp_Closest'))
                resamp_sum_interp = np.sum(grid.multiply('Land Fraction','Resamp_Interp'))

                plt.rcParams.update({'font.size': 18})
                from matplotlib import cm
                cmap = cm.winter

                diff_ideal = land_water_contrast*(target_sum-resamp_sum)
                diff_closest = land_water_contrast*(target_sum_shift-resamp_sum_closest)
                diff_interp = land_water_contrast*(target_sum_shift-resamp_sum_interp)

                diff_list_ideal.append(diff_ideal)   
                diff_list_closest.append(diff_closest)   
                diff_list_interp.append(diff_interp)
                distance_list_closest.append(distance_closest)   
                num_done += 1
                
                out_str = f'{scene}:{band_name[band_num+freq_index_offset]}: num = {num_done} isamp= {isamp} lon = {lon:.2f}  lat = {lat:.2f}'
                out_str = f'{out_str} Ideal: {diff_ideal:.3f}'
                out_str = f'{out_str} Closest: {diff_closest:.3f}'
                out_str = f'{out_str} Interp: {diff_interp:.3f}'
                out_str = f'{out_str} Distance: {distance_closest:.3f}'
                print(out_str)
            
            diff_list_ideal = np.array(diff_list_ideal)
            diff_list_closest = np.array(diff_list_closest)
            diff_list_interp = np.array(diff_list_interp)
            distance_closest = np.array(distance_list_closest)

            diff_list_ideal_all[band_num,scene_index,:] = diff_list_ideal
            diff_list_closest_all[band_num,scene_index,:] = diff_list_closest
            diff_list_interp_all[band_num,scene_index,:] = diff_list_interp
            distance_closest_all[band_num,scene_index,:] = distance_list_closest

            print('------------------------------------------------------------------')
            print(f'{scene}:{band_name[band_num+freq_index_offset]} Beamwidth={beamwidth}')
            print(f'Ideal:   {np.mean(diff_list_ideal):.4f} +/- {np.std(diff_list_ideal):.4f} ')
            print(f'Closest:  {np.mean(diff_list_closest):.4f} +/- {np.std(diff_list_closest):.4f} ')
            print(f'Interp:  {np.mean(diff_list_interp):.4f} +/- {np.std(diff_list_interp):.4f} ')
            print('------------------------------------------------------------------')
            lons_all = np.array(lons_all)
            lats_all = np.array(lats_all)

            # if beamwidth == 60:
            #     txt_file = plot_path / f'random_locations_{scene}_{beamwidth}km.txt'
            #     with open(txt_file,'w') as f:
            #         for i in range(len(lons_all)):
            #             if scene in ['lakes','midwest','coastline']:
            #                 f.write(f'{lons_all[i]:.4f}, {lats_all[i]:.4f}, {shift_x_all[i]:.4f},  {shift_y_all[i]:.4f}\n')
            #             else:
            #                 f.write(f'{lons_all[i]:.4f}, {lats_all[i]:.4f}, {shift_x_all[i]:.4f},  {shift_y_all[i]:.4f}, {phi_all[i]:.4f}\n')


            lons_all[lons_all > 180.0] = lons_all[lons_all > 180.0] - 360        
            
            std_diff_table_ideal[band_num,scene_index] = np.std(diff_list_ideal)
            std_diff_table_closest[band_num,scene_index] = np.std(diff_list_closest)
            std_diff_table_interp[band_num,scene_index] = np.std(diff_list_interp)
            
            rms_diff_table[band_num,scene_index,0] = rms(diff_list_ideal)
            rms_diff_table[band_num,scene_index,1] = rms(diff_list_interp)
            rms_diff_table[band_num,scene_index,2] = rms(diff_list_closest)

            bin_edges = [0.0,0.5,1.0,1.5,2.0,2.5,3.5]

            binned_rms = np.zeros((6))
            num_in_bin = np.zeros((6))

            # for i in range(0,6):
            #     diffs_in_bin_closest = diff_list_closest[np.all([distance_list_closest > bin_edges[i],distance_list_closest <= bin_edges[i+1]],axis=0)]
            #     diffs_in_bin_interp = diff_list_interp[np.all([distance_list_closest > bin_edges[i],distance_list_closest <= bin_edges[i+1]],axis=0)]
            #     print(f'{bin_edges[i]:.2f}-{bin_edges[i+1]:.2f}   {rms(diffs_in_bin_closest):.2f}  {rms(diffs_in_bin_interp):.2f} {len(diffs_in_bin_closest)}')
            

            fig,ax = plt.subplots()

    rms_diff_table_xr = xr.Dataset(data_vars = dict(rms_diff = (["bands","scenes","colloc_types"],rms_diff_table),
                                                    diff_ideal_all = (["bands","scenes","samples"],diff_list_ideal_all),
                                                    diff_closest_all = (["bands","scenes","samples"],diff_list_closest_all),
                                                    diff_interp_all = (["bands","scenes","samples"],diff_list_interp_all),
                                                    distance_closest_all = (["bands","scenes","samples"],distance_closest_all)),

                                   coords = dict(bands = (["bands"],np.arange(len(band_list))),
                                                 band_names = (["band_names"],[band_name[i+freq_index_offset] for i in band_list]),

                                                 scenes = (["scenes"],np.arange(len(scene_list))),
                                                 scene_names = (["scene_names"],scene_list),

                                                 colloc_types = (["colloc_types"],np.arange(3)),
                                                 colloc_names = (["colloc_names"],["ideal","interp","closest"]),
                                                 
                                                 samples = (["samples"],np.arange(num_to_do)),
                                                 ),
                                   attrs=dict(sensor_name = sensor_name,
                                              ksat = ksat,
                                              beamwidth = beamwidth,
                                              land_water_contrast = land_water_contrast,
                                              gradient_magnitude =gradient_magnitude,  #K/kilometer or 50K / 100 km
                                              kscan = kscan,
                                              fov = fov,
                                              num_to_do = num_to_do)
                                    )
    nc_file = plot_path / f'rms_diff_table__and_differnce_lists_{sensor_name}_{ksat}_{beamwidth}km_fov_{fov}_kscan_{kscan}_xr.nc'
    rms_diff_table_xr.to_netcdf(nc_file)
    print()
    print('---------------------------------------------------')
    print(f'sat_name = {sensor_name}')
    if sensor_name in ['SSMIS']:
        print(f'ksat = {ksat}')
    print(f'land_water_contrast = {land_water_contrast}')
    print(f'kscan = {kscan}')
    print(f'fov = {fov}')
    print(f'beamwidth = {beamwidth}')
    print()
    print('RMS Difference Ideal')
    print(scene_list)
    for band_index,band_num in enumerate(band_list):
        out_str = f'{band_name[band_num+freq_index_offset]}'
        for scene_index,scene in enumerate(scene_list):
            out_str = f'{out_str}   {rms_diff_table[band_index,scene_index,0]:.5f}'
        print(out_str)

    print()
    print('RMS Difference Interpolated')
    print(scene_list)
    for band_index,band_num in enumerate(band_list):
        out_str = f'{band_name[band_num+freq_index_offset]}'
        for scene_index,scene in enumerate(scene_list):
            out_str = f'{out_str}   {rms_diff_table[band_index,scene_index,1]:.5f}'
        print(out_str)


    print()
    print('RMS Difference closest')
    print(scene_list)
    for band_index,band_num in enumerate(band_list):
        out_str = f'{band_name[band_num+freq_index_offset]}'
        for scene_index,scene in enumerate(scene_list):
            out_str = f'{out_str}   {rms_diff_table[band_index,scene_index,2]:.5f}'
        print(out_str)

    # EXTRA code, previously used for plotting


            
        # grid.addlayer(name='Target_Shift')
        # grid.setdata(target_sub_shift,name='Target_Shift')
        # grid.normalize(name='Target_Shift')

        # grid.addlayer(name='Resample')
        # grid.setdata(resample_sub_array[:,:,0,0],name='Resample')
        # grid.normalize(name='Resample')

        # grid.addlayer(name='Resample_Interp')
        # grid.setdata(resample_sub_interp,name='Resample_Interp')
        # grid.normalize(name='Resample_Interp')

        # diff = grid.subtract('Target','Resample')
        # grid.addlayer(name = 'Diff')
        # grid.setdata(diff,'Diff')

        # diff_shifted = grid.subtract('Target_Shift','Resample')
        # grid.addlayer(name = 'Diff_Shifted')
        # grid.setdata(diff_shifted,'Diff_Shifted')

        # diff_shifted_interp = grid.subtract('Target_Shift','Resample_Interp')
        # grid.addlayer(name = 'Diff_Shifted_Interp')
        # grid.setdata(diff_shifted_interp,'Diff_Shifted_Interp')


        # maxval = np.max(target_sub)
        # cbticks = np.array([-0.001,-0.0005,0.0,0.0005,0.001])
        # fig,ax,im = grid.contourplot(name='Target',title ='30 km Target',vmin=-maxval,vmax=maxval,scale=False,cbticks=cbticks)
        # png_file = plot_path / 'target.png'
        # fig.savefig(png_file)

        # fig,ax,im = grid.contourplot(name='Target_Shift',title ='30 km Target, shifted diagonally ~2.8 km',vmin=-maxval,vmax=maxval,scale=False,cbticks=cbticks)
        # png_file = plot_path / 'target_shift.png'
        # fig.savefig(png_file)
        
        # fig,ax,im = grid.contourplot(name='Resample',title ='Closest Resampled',vmin=-maxval,vmax=maxval,scale=False,cbticks=cbticks)
        # png_file = plot_path / 'resamp_closest.png'
        # fig.savefig(png_file)
        
        # fig,ax,im = grid.contourplot(name='Resample_Interp',title ='Resampled, Interpolated',vmin=-maxval,vmax=maxval,scale=False,cbticks=cbticks)
        # png_file = plot_path / 'resamp_interp.png'
        # fig.savefig(png_file)
        
        # fig,ax,im = grid.contourplot(name='Diff',title ='Pattern Difference, Ideal',vmin=-maxval/10.,vmax=maxval/10.,scale=False,cbticks=cbticks/10.)
        # png_file = plot_path / 'diff_ideal.png'
        # fig.savefig(png_file)
        
        # fig,ax,im = grid.contourplot(name='Diff_Shifted',title ='Pattern Difference, Closest',vmin=-maxval/10.,vmax=maxval/10.,scale=False,cbticks=cbticks/10.)
        # png_file = plot_path / 'diff_closest.png'
        # fig.savefig(png_file)
        
        # fig,ax,im = grid.contourplot(name='Diff_Shifted_Interp',title ='Pattern Difference, Interpolated',vmin=-maxval/10.,vmax=maxval/10.,scale=False,cbticks=cbticks/10.)
        # png_file = plot_path / 'diff_interp.png'
        # fig.savefig(png_file)
