# I used this to determine the scale factor in target_gain in 
# antenna\AMSRE\AMSR2_Antenna_Gain.py.  The scale factor sets the 
# size of the footprint in km space.  CM July 23 2021

import numpy as np 
from AMSR2_Antenna_Gain import *
import matplotlib.pyplot as plt


dist_km = np.arange(0,20000)*0.01

gains = target_gain(dist_km,diameter_in_km = 30.0)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.plot(dist_km,gains)

ok = gains < 0.5
print(dist_km[ok][0])
print(dist_km[ok][0]/15.0)
plt.show()
print
