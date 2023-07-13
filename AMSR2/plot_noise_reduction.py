
import numpy as np 
import xarray as xr 
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import cm



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

noise_reduction_factor = np.zeros((485,4))

for band in [1,2,3]:
    weights = read_footprint_weight_file(
                sensor = 'AMSR2',
                target_diameter = 30,
                band = band,
                data_root = 'L:/access/resampling/',
                verbose = True)


    for fov in np.arange(0,485):
        slice  = weights[1,fov,:,:]
        noise_reduction_factor[fov,band] = np.sqrt(np.sum(slice*slice))


fig,ax = plt.subplots(figsize=(10,6))
fovs = np.arange(1,486)
ax.plot(fovs,noise_reduction_factor[:,1],label='11 GHz')
ax.plot(fovs,noise_reduction_factor[:,2],label='19 GHz')
ax.plot(fovs,noise_reduction_factor[:,2],label='24 GHz')
ax.set_ylim(0.0,1.0)
ax.set_xlabel('Resampled Field of View',fontsize=14)
ax.set_ylabel('Noise Reduction Factor',fontsize=14)
ax.legend()
plt.show()

print