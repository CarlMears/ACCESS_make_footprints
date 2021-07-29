#     From O:\amsr_l2\v05\resampling\show_antenna_patterns.f


#     The values of these AMSR-E parameters come from E:\AMSR\Resampling_weights\PM1_AMSR_MOD.for
#     The values of th Midori2   parameters come from E:\AMSR\Resampling_weights\ADEOSII_AMSR_MOD.for
#     For these antenna parameters, Peter found one set for 89 GHZ, which was the averaged of 89A and 89B
#     This value was put in the IFREQ=8 slot
#     The IFREQ=6 and 7 slots were set to 0
#     For Midori2, the IFREQ=6 and 7 slots contain the 50.3 and  52.8 patterns respectively.

#     Here, I use slot 7 for 89A and slot 8 for 89B, so they have the same value
#     When I do Midori-2, I will need to decide what to do about the 50 GHz channels (maybe put in average value into slot 6)
  
#  converted to python so I can see what's going one

import numpy as np
import matplotlib.pyplot as plt
from AMSR2_Antenna_Gain import *


sin_anglim   = np.array([0.105,    0.105,     0.0524,   0.0524,   0.0262,   0.000,  0.0262 ,    0.0262],dtype=np.float64)
ant_approx_a = np.array([4.343e-6, 2.096e-6,	1.890e-6, 1.623e-6, 7.248e-7, 0.000,  2.070e-6,   2.070e-6],dtype=np.float64)
ant_approx_b = np.array([6.892e-4, 4.059e-4,	3.727e-4, 7.251e-4, 3.051e-4, 0.000,  2.381e-4,   2.381e-4],dtype=np.float64)
ant_approx_c = np.array([0.503,    0.662,     1.391,    1.804,    1.964,    0.000,  4.593,      4.593],dtype=np.float64)
ant_approx_d = np.array([0.651,    1.345, 	4.844,    4.721,    15.18,    0.000,  79.785,     79.785],dtype=np.float64)

jaxa_width = [1.8,1.2,0.65,0.75,0.35]
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

for jfreq in range(0,5):
    anglim =   np.rad2deg(np.arcsin(sin_anglim[jfreq]))
    coeff_a =  ant_approx_a[jfreq]
    coeff_b =  ant_approx_b[jfreq]
    coeff_c =  ant_approx_c[jfreq]
    coeff_d =  ant_approx_d[jfreq]

    #print(f"anglim = {anglim:.3f}")

    delta = np.arange(0,70001).astype(np.float64)*0.0001
    sind_delta = np.sin(np.deg2rad(delta))

    #gain = coeff_a + coeff_b*np.exp(-coeff_c*delta) + 	np.exp(-coeff_d*delta*delta)
    gain = AMSRE_antenna_gain(delta,jfreq)
    mx = gain[0]
    half_max_angle = np.nanmin(delta[gain < 0.5*mx])
    print(jfreq,2.0*half_max_angle,1.0/(2.0*half_max_angle/jaxa_width[jfreq]))

    ax.plot(delta,gain)
    #ax.plot(delta,gain2,color='blue')
    #ax.plot(delta,gain*sind_delta,color='red')

    gain[delta > anglim] = 1.0e-6

    qsum = np.sum(gain*sind_delta)
    sumgain = qsum*2.0*np.pi*np.deg2rad(0.0001)
ax.set_xlim([0.0,3.0])
plt.show()
#         ok = delta <= anglim

# 	    QSUM=0; DELTASV=-1
# 	    DO I=0,70000
# 	        DELTA=0.0001D0*I

#         IF (DELTA.LE.ANGLIM) THEN
# 	        gain = coeff_a + coeff_b*exp(-coeff_c*delta) + 	exp(-coeff_d*delta*delta)
# 	        IF(I.EQ.0) 
#               GAIN0=GAIN
# 	        else
# 	            gain=1.e-6
# 	        endif

# 	    WRITE(4) GAIN

# 	QSUM=QSUM + GAIN*SIND(DELTA)
# 	IF(ABS(GAIN/GAIN0 -0.5).LT.0.001) DELTASV=DELTA
# 	ENDDO !I

# 	SUMGAIN=QSUM*6.283185D0*1.745329D-6	  !TWO PI TIMES INTEGRATION STEP IN RADIANS (0.001 DEG)

# 	WRITE(6,1001) JFREQ, SUMGAIN, 2*DELTASV,
#      & 10*DLOG10((coeff_a + coeff_b*exp(-coeff_c*ANGLIM) + exp(-coeff_d*ANGLIM*ANGLIM))/GAIN0) 
#  1001 FORMAT(I3,E15.5,F8.3,F8.2)
# 	ENDDO !JFREQ


	 					 		 		     
