def read_gain_array(sat_name='AMSR2',beamwidth=30,band_num=4,kscan=1,fov=1,gain_type='source',verbose = False):
    import numpy as np
    import struct

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

    elif gain_type == 'resample_v3':
        path = f"L:access/resampling/{sat_name}/resample_gains_v3/circular_{beamwidth:02d}km/band_{band_num:02d}/"
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

                iy,ix,gain = struct.unpack('<hhd',x)
              
                #if num_lines_read < 3:
                #    print(iy,ix,gain)
                num_lines_read += 1
                gain_array[iy,ix]  =gain
                gain_tot += gain
                if gain > max_gain:
                    max_gain = gain
                    iy_max = iy
                    ix_max = ix
                
            except:
                continue
    if verbose:
        print(f'Num Lines read = {num_lines_read}')
        print(f'Total Gain: {gain_tot:8.3f}')
        print()
        print(iy_max,ix_max)

    return gain_array,iy_max,ix_max
