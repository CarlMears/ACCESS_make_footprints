# read in footprint weight file

import numpy as np 
import xarray as xr 
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import cm

from rss_plotting.plot_2d_array import plot_2d_array

def read_footprint_weight_file(
                sensor = 'AMSR2',
                target_diameter = 30,
                band = 1,
                data_root = 'L:/access/resampling/',
                verbose = False):

    data_root = Path(data_root) 
    file = data_root / f'{sensor}/resample_weights_v3/circular_{target_diameter:02d}km/resample_weights_band_{band:02d}.dat'

    if verbose:
        print(f'Reading: {file}')

    with open(file,mode='rb') as f:
        if sensor == 'AMSR2':
            dims = (2,485,29,243)
        else:
            raise ValueError(f'Sensor {sensor} not implemented')
            
        tot_size = dims[0]*dims[1]*dims[2]*dims[3]
        weights = np.fromfile(file, dtype=np.float64, count=tot_size).reshape(dims)

    return weights


if __name__ == '__main__':

    weights = read_footprint_weight_file()
    print(weights.shape)

    fov = 20
    scan = 0

    slice = weights[scan,fov,:,:]

    print(np.sum(slice))
    print(np.sum(slice*slice))
    print(np.sqrt(np.sum(slice*slice)))
    print
