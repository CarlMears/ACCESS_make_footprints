def _read_gain_array_AMSR2(beamwidth=30,band_num=4,kscan=1,fov=1,gain_type='source',verbose = False):
    import numpy as np
    import struct
    from pathlib import Path

    gain_array = np.zeros((1602,2402),dtype = 'float64')
    sat_name = 'AMSR2'

    if gain_type == 'target':
        path = f"L:/access/resampling/{sat_name}/target_gains/{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'target_v2':
        path = f"L:/access/resampling/{sat_name}/target_gains_v2/{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'source':
        path = f"L:/access/resampling/{sat_name}/source_gainsband_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'source_v2':
        path = f"L:/access/resampling/{sat_name}/source_gains_v2/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'source_v3':
        path = f"L:/access/resampling/{sat_name}/source_gains_v3/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:03d}.dat"

    elif gain_type == 'resample':
        path = f"L:/access/resampling/{sat_name}/resample_gains/circular_{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:04d}.dat"

    elif gain_type == 'resample_v2':
        path = f"L:/access/resampling/{sat_name}/resample_gains_v2/circular_{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:04d}.dat"

    elif gain_type == 'resample_v3':
        path = f"L:/access/resampling/{sat_name}/resample_gains_v3/circular_{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = f"{path}s{kscan:02d}c{fov:04d}.dat"
    else:
        raise ValueError(f'gain type {gain_type} makes no sense')

    num_lines_read = 0
    gain_tot = 0.0
    max_gain = 0.0
    path_to_gain_file = Path(file_name)
    num_pts=int(path_to_gain_file.stat().st_size/12)
    if verbose:
        print(f'Reading {num_pts} points from {path_to_gain_file}')
    with open(path_to_gain_file,'rb') as f:
        for i in range(0,num_pts):
            try:
                x = f.read(12)

                iy,ix,gain = struct.unpack('<hhd',x)
              
                # if gain_type == 'resample_v3':
                #     print(iy,ix,gain)
                #     pass

                #if num_lines_read < 3:
                #    print(iy,ix,gain)
                num_lines_read += 1
                if np.isfinite(gain):
                    gain_array[iy,ix]  =gain
                    gain_tot += gain
                    if gain > max_gain:
                        max_gain = gain
                        iy_max = iy
                        ix_max = ix
                else:
                    print('not Finite')
                
            except:
                continue
    if verbose:
        print(f'Num Lines read = {num_lines_read}')
        print(f'Total Gain: {gain_tot:8.3f}')
        print()
        print(iy_max,ix_max)

    return gain_array,iy_max,ix_max

def _read_gain_array_SSMIS(ksat=18,beamwidth=30,band_num=4,kscan=1,fov=1,gain_type='source',verbose = False):
    import numpy as np
    import struct
    from pathlib import Path

    gain_array = np.zeros((1602,2402),dtype = 'float64')
    
    root_path = Path(f"L:/access/resampling/SSMIS/f{ksat:02d}")
    if gain_type == 'target':
        path = root_path / f"target_gains/{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = path / f"s{kscan:02d}c{fov:03d}.dat"
    elif gain_type == 'source':
        path = root_path / f"source_gains/band_{band_num:02d}/"
        file_name = path / f"s{kscan:02d}c{fov:03d}.dat"
    elif gain_type == 'resample':
        path = root_path / f"resample_gains/{beamwidth:02d}km/band_{band_num:02d}/"
        file_name = path / f"s{kscan:02d}c{fov:03d}.dat"
    else:
        raise ValueError(f'gain type {gain_type} makes no sense')

    num_lines_read = 0
    gain_tot = 0.0
    max_gain = 0.0

    path_to_gain_file = Path(file_name)
    num_pts=int(path_to_gain_file.stat().st_size/12)
    if verbose:
        print(f'Reading {num_pts} points from {path_to_gain_file}')
    with open(path_to_gain_file,'rb') as f:
        for i in range(0,num_pts):
            try:
                x = f.read(12)
                iy,ix,gain = struct.unpack('<hhd',x)
                num_lines_read += 1
                if np.isfinite(gain):
                    gain_array[iy,ix]  =gain
                    gain_tot += gain
                    if gain > max_gain:
                        max_gain = gain
                        iy_max = iy
                        ix_max = ix
                else:
                    print('not Finite')
 
            except:
                continue
    if verbose:
        print(f'Num Lines read = {num_lines_read}')
        print(f'Total Gain: {gain_tot:8.3f}')
        print()
        print(iy_max,ix_max)

    return gain_array,iy_max,ix_max

def _read_gain_array_SSMI(ksat=13,beamwidth=30,freq_num=0,kscan=1,fov=1,gain_type='source',verbose = False):
    import numpy as np
    import struct
    from pathlib import Path

    gain_array = np.zeros((1602,2402),dtype = 'float64')
    
    root_path = Path(f"L:/access/resampling/SSMI/f{ksat:02d}")
    if gain_type == 'target':
        path = root_path / f"target_gains/{beamwidth:02d}km/freq_{freq_num:02d}/"
        file_name = path / f"s{kscan:02d}c{fov:03d}.dat"
    elif gain_type == 'source':
        path = root_path / f"source_gains/freq_{freq_num:02d}/"
        file_name = path / f"s{kscan:02d}c{fov:03d}.dat"
    elif gain_type == 'resample':
        path = root_path / f"resample_gains/{beamwidth:02d}km/freq_{freq_num:02d}/"
        file_name = path / f"s{kscan:02d}c{fov:03d}.dat"
    else:
        raise ValueError(f'gain type {gain_type} makes no sense')

    num_lines_read = 0
    gain_tot = 0.0
    max_gain = 0.0

    path_to_gain_file = Path(file_name)
    num_pts=int(path_to_gain_file.stat().st_size/12)
    if verbose:
        print(f'Reading {num_pts} points from {path_to_gain_file}')
    with open(path_to_gain_file,'rb') as f:
        for i in range(0,num_pts):
            try:
                x = f.read(12)
                iy,ix,gain = struct.unpack('<hhd',x)
                num_lines_read += 1
                if np.isfinite(gain):
                    gain_array[iy,ix]  =gain
                    gain_tot += gain
                    if gain > max_gain:
                        max_gain = gain
                        iy_max = iy
                        ix_max = ix
                else:
                    print('not Finite')
 
            except:
                continue
    if verbose:
        print(f'Num Lines read = {num_lines_read}')
        print(f'Total Gain: {gain_tot:8.3f}')
        print()
        print(iy_max,ix_max)

    return gain_array,iy_max,ix_max

def read_gain_array(*,sat_name,beamwidth,band_num,kscan,fov,gain_type='source',verbose = False,ksat=18):

    if sat_name == 'AMSR2':
        if gain_type == 'source':
            gain_type = 'source_v3'
        if gain_type == 'target':
            gain_type = 'target_v2'
        gain_array,iy_max,ix_max = _read_gain_array_AMSR2(beamwidth=beamwidth,
                                                          band_num=band_num,
                                                          kscan=kscan,
                                                          fov=fov,
                                                          gain_type=gain_type,
                                                          verbose = verbose)
    elif sat_name == 'SSMIS':
        gain_array,iy_max,ix_max = _read_gain_array_SSMIS(ksat=ksat,
                                                          beamwidth=beamwidth,
                                                          band_num=band_num,
                                                          kscan=kscan,
                                                          fov=fov,
                                                          gain_type=gain_type,
                                                          verbose = verbose)
    elif sat_name == 'SSMI':
        gain_array,iy_max,ix_max = _read_gain_array_SSMI(ksat=ksat,
                                                          beamwidth=beamwidth,
                                                          freq_num=band_num,
                                                          kscan=kscan,
                                                          fov=fov,
                                                          gain_type=gain_type,
                                                          verbose = verbose)
    else:
        raise ValueError(f'sat_name = {sat_name} invalid')

    return gain_array,iy_max,ix_max

    

