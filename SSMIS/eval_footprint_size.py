import numpy as np
import matplotlib.pyplot as plt

from read_gain_array import read_gain_array_SSMIS

kscan = 28
fov = 25
band_num=3
ksat = 18

gain_array,iy_max,ix_max = read_gain_array_SSMIS(ksat=ksat,
                                           band_num=band_num,
                                           kscan=kscan,
                                           fov=fov,
                                           gain_type='source',
                                           verbose = True)
size1 = np.sum(gain_array[iy_max-100:iy_max+100,ix_max] > 0.5)
size2 = np.sum(gain_array[iy_max,ix_max-100:ix_max+100] > 0.5)
print(f'Approx Source Footprint Size (band={band_num}) = {size1:.1f} x {size2:.1f}')

gain_array,iy_max,ix_max = read_gain_array_SSMIS(ksat=ksat,
                                           beamwidth=70,
                                           band_num=2,
                                           kscan=1,
                                           fov=200,
                                           gain_type='target',
                                           verbose = True)
size1 = np.sum(gain_array[iy_max-100:iy_max+100,ix_max] > 0.5)
size2 = np.sum(gain_array[iy_max,ix_max-100:ix_max+100] > 0.5)
print(f'Approx Target Footprint Size (band={band_num}) = {size1:.1f} x {size2:.1f}')                                 
print







