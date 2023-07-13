import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
from matplotlib import cm

#local packages specific to this work
from read_SSMI_native_location_file import read_SSMI_native_location_file
from SSMI_Antenna_Gain import SSMI_antenna_gain
from rss_gridding.local_earth_grid import LocalEarthGrid

def calc_delta_theta_y(range,incidence,dy):
    import numpy as np 
    inc_rad = np.deg2rad(incidence)
    r1 = np.sqrt(range*range + 
                 dy*dy + 
                 2*range*dy*np.sin(np.deg2rad(incidence)))

    z = (np.square(dy) - np.square(r1-range))/(4*range*r1)
    return 2.0*np.rad2deg(np.arcsin(np.sqrt(z))),r1

def calc_delta_theta_x(slant_range,dx):
    import numpy as np
    return np.rad2deg(np.arcsin(dx/slant_range))

def calc_delta_theta(slant_range,theta,azim,x,y,max_distance_to_consider):

    delta_theta = np.ones_like(x)*90.0
    x1 =  np.cos(np.deg2rad(azim))*x - np.sin(np.deg2rad(azim))*y
    y1 =  np.sin(np.deg2rad(azim))*x + np.cos(np.deg2rad(azim))*y

    delta_theta_y,r1 = calc_delta_theta_y(slant_range,theta,y1)
    delta_theta_x = np.rad2deg(np.arcsin(x1/slant_range))

    ok = np.sqrt(np.square(x1)+np.square(y1)) < max_distance_to_consider
    delta_theta[ok] = (np.sqrt(np.square(delta_theta_y) + np.square(delta_theta_x)))[ok]
    return delta_theta,r1/slant_range

def write_footprint_locations(lats,lons,azim,ksat,footprint_type,freq = 2):
    import os

    if footprint_type == 'target':
        textfile = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/target_gains/locs/find_target_gain_f{ksat:02d}.loc')
    elif footprint_type == 'source':
        textfile = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/source_gains/locs/find_source_gain_f{ksat:02d}_freq_{freq:02d}_locs.txt')
    else:
        raise ValueError(f'Invalid Source Type: {footprint_type}')
    
    sz = lats.shape
    os.makedirs(textfile.parent,exist_ok=True)

    with open(textfile,'w') as f:
        for iscan in range(0,sz[0]):
            for ifov in range(0,sz[1]):
                if footprint_type == 'target':
                    f.write(f' {iscan:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f}\n')
                else:
                    print(f' {freq:04d} {iscan+1:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f} {azim[iscan,ifov]:12.5f}\n')
                    f.write(f' {freq:04d} {iscan+1:04d} {ifov+1:04d} {lats[iscan,ifov]:12.5f} {lons[iscan,ifov]:12.5f} {azim[iscan,ifov]:12.5f}\n')
                    



    
ksat = 13
orbit_num = 5000
scan0 = 456
scan0_hi = 912

freq_list = ['19','22','37','85']
make_plots = True
footprint_type = 'source'

locs = read_SSMI_native_location_file(ksat,orbit_num)

# get the subset of lats and lons
for ifreq,freq in enumerate(freq_list):
    if ifreq < 3:
        delta_scans = 14
        scan0 = 456
    else:
        delta_scans = 28
        scan0 = 912
    
    output_path_binary = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/source_gains_test/')
    output_path_plot = Path(f'L:/access/resampling/SSMI/f{ksat:02d}/source_gains_test/plots/')
    
    os.makedirs(output_path_binary,exist_ok=True)
    os.makedirs(output_path_plot,exist_ok=True)

    #find the mean values for slant_range and incidence angle
    rng = locs[f'rng{freq}']
    tht = locs[f'tht{freq}']

    slant_range = np.mean(rng[rng > 10000.0])
    theta = np.mean(tht[tht > 10.0])

    lons = locs[f'lons{freq}']
    lats = locs[f'lats{freq}']
    azim = locs[f'phi{freq}']

    num_cells = lons.shape[1]

    lats = lats[scan0-delta_scans:scan0+delta_scans+1,:]
    lons = lons[scan0-delta_scans:scan0+delta_scans+1,:]-lons[scan0,int(num_cells/2)]
    azim = azim[scan0-delta_scans:scan0+delta_scans+1,:]

    write_footprint_locations(lats,lons,azim,ksat,footprint_type,freq = ifreq)

    #1km spacing very close to Frank's 0.01 degree spacing
    grid = LocalEarthGrid(center_lon = 0.0,center_lat = 0.0,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 2401,num_y = 1601,
                          name='Single')

    if make_plots:
        grid.addlayer(name='Sum')

    x = grid.xlocs*1000.0
    y = grid.ylocs*1000.0

    num_scans = lats.shape[0]
    num_fovs = lats.shape[1]

    for jscan in np.arange(0,num_scans):
        for jfov in np.arange(0,num_fovs):
        #for jfov in np.arange(0,34):
            xloc,yloc = grid.projection(lons[jscan,jfov],lats[jscan,jfov])
            print(freq,jscan-delta_scans,jfov,lons[jscan,jfov],lats[jscan,jfov],xloc,yloc)

            # calculate the cross-scan vector in the x,y plane in lat/lon units.
            if (jfov == 0):
                cross_scan_vector = 0.5*np.array([lons[jscan,1] - lons[jscan,0],
                                                    lats[jscan,1] - lats[jscan,0]])
            elif (jfov == num_fovs-1):
                cross_scan_vector = 0.5*np.array([lons[jscan,num_fovs-1] - lons[jscan,num_fovs-2],
                                                    lats[jscan,num_fovs-1] - lats[jscan,num_fovs-2]])
            else:
                cross_scan_vector = 0.25*np.array([lons[jscan,jfov+1] - lons[jscan,jfov-1],
                                                    lats[jscan,jfov+1] - lats[jscan,jfov-1]])

            gain_tot = np.zeros((1601,2401))
            # This loop smears the footprint in the cross-track direction to account for the
            # non-zero sampling time of the SSMI instrument
            for smear_factor in [-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9]:
                xloc,yloc = grid.projection(lons[jscan,jfov]+smear_factor*cross_scan_vector[0],
                                            lats[jscan,jfov]+smear_factor*cross_scan_vector[1])

                # calculate delta theta from boresight vector.
                max_distance_to_consider = 200000.0
                delta_theta,r_norm = calc_delta_theta(slant_range,theta,
                                            azim[jscan,jfov],
                                            x-xloc,y-yloc,
                                            max_distance_to_consider)

                # calculate 2-d antenna gain 
                gain = SSMI_antenna_gain(delta_theta,ifreq)/np.square(r_norm)
                gain_tot = gain_tot + gain
            gain = gain_tot/7.0
            
            # store it in the local grid instance
            grid.setdata(gain,name='Single')
 
            if make_plots:
                if ((jfov%8 == 4) and (jscan%(delta_scans/2) == 0)):
                    print(f'Adding scan {jscan} fov {jfov} to plot')
                    grid.adddata(gain,name='Sum')
            
            output_path_binary2 = output_path_binary / f'freq_{ifreq:02d}/'
            os.makedirs(output_path_binary2,exist_ok = True)
            source_file_name = output_path_binary2 /f's{jscan+1:02d}c{jfov+1:03d}.dat'
            grid.write_binary(outfile = source_file_name,name='Single',save_thres=0.0001) 
            # note that this method only stores the data if it is above the threshold value
# plot it if desired
    if make_plots:
        fig0 = plt.figure(figsize = (10,11))
        outfile = output_path_plot / f'decimated_scan_{freq}.png'
        os.makedirs(outfile.parent,exist_ok = True)
        fig0,ax,im =grid.contourplot(name='Sum',
                        units='Footprint Weight',
                        title=f'{footprint_type} footprints, {freq_list[ifreq]}',
                        fig_in=fig0,
                        vmin=0.0,vmax=1.0,
                        cmap = cm.viridis,
                        plt_contours = False,
                        scale=False,
                        cbticks = np.array([0.0,0.2,0.4,0.6,0.8,1.0]),
                        outfile = outfile)
        plt.show()

    print('Normal End')



