import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


plt.rcParams.update({'font.size': 16})
sensor_name = 'SSMIS'

if sensor_name == 'SSMIS':
    #SSMIS
    num_target_scans = 2
    num_target_fovs = 446
    num_source_scans = 57
    num_source_fovs = 90

    path_to_weights = Path('L:/access/resampling/SSMIS/f18/resample_weights/70km')
    plt_path = Path('L:/access/resampling/SSMIS/f18/resample_gains/70km/plots')
    band_name = ['','','19GHz','22GHz','37GHz']
elif sensor_name == 'AMSR2':
    
    num_target_scans = 2
    num_target_fovs = 446
    num_source_scans = 57
    num_source_fovs = 90

    path_to_weights = Path('L:/access/resampling/AMSR2/resample_weights/70km')
    plt_path = Path('L:/access/resampling/SSMIS/f18/resample_gains/70km/plots')
    band_name = ['7GHz','11GHz','19GHz','22GHz','37GHz']
else:
    raise ValueError(f'sensor_name = {sensor_name} not valid')




fig,ax = plt.subplots(figsize=(10,6))
for band in [2,4]:

    weight_file = path_to_weights / f'resample_weights_band_{band:02d}.dat'

    #statement in fortran
    #real(real64) :: weights(num_source_fovs,num_source_scans,num_target_fovs,num_target_scans)

    total_size_should_be = num_target_scans*num_target_fovs*num_source_scans*num_source_fovs
    print(total_size_should_be)

    weights = np.fromfile(weight_file)
    print(len(weights))

    weights = np.reshape(weights,(num_target_scans,num_target_fovs,num_source_scans,num_source_fovs))

    noise_factor = np.zeros((num_target_scans,num_target_fovs))

    for iscan in range(num_target_scans):
        for ifov in range(num_target_fovs):
            noise_factor[iscan,ifov] = np.sqrt(np.sum(np.square(weights[iscan,ifov,:,:])))



    ax.plot(noise_factor[0,:],label = band_name[band])


ax.set_ylim(0.0,0.4)
ax.legend()
ax.set_xlabel('FOV Number')
ax.set_ylabel('Noise Reduction Factor')
ax.set_title('SSMIS')
ax.tick_params(direction = 'in')

png_file = plt_path / f'resampled_noise_factor.png'
fig.savefig(png_file)
plt.show()
print