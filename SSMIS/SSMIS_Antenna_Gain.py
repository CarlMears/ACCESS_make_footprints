#     From O:\ssmis\resampling\find_gains.f

#     These are Frank's comments on the fortran code.
#     !see memo2 in O:\ssmis\resampling

#  Comments below by CM, Summer 2022
#
#  converted to python so I can see what's going on
#  the python version is used to create the input data for 
#  the resampling weight calculations
#

 
def SSMIS_antenna_gain(delta,band):
    import numpy as np

    theta_full   = np.array([0.00,      0.00,     2.00,     2.00,     1.22,    0.00,     0.00])
    theta_half   = theta_full/2.0
    sin_anglim   = np.array([0.00,      0.00,     0.0524,   0.0524,   0.0262,   0.000,    0.000],dtype=np.float64)
    ant_approx_a = np.array([0.00,      0.00,     2.00e-6,  2.00e-6,  2.00e-6,  2.00e-6,  2.00e-6],dtype=np.float64)
    ant_approx_b = np.array([0.00,      0.00,     4.00e-4,  4.00e-4,  4.00e-4,  4.00e-4,  4.00e-4],dtype=np.float64)
    ant_approx_c = np.zeros(6,dtype=np.float64)
    ant_approx_d = np.zeros(6,dtype=np.float64)
    ant_approx_d[2] = -np.log(0.5)/theta_half[2]**2
    ant_approx_d[3] = -np.log(0.5)/theta_half[3]**2
    ant_approx_d[4] = -np.log(0.5)/theta_half[4]**2

    ant_approx_c[2] = np.sqrt(0.3*ant_approx_d[2])
    ant_approx_c[3] = np.sqrt(0.3*ant_approx_d[3])
    ant_approx_c[4] = np.sqrt(0.3*ant_approx_d[4])
        #                          0         1          2         3         4        5       6            7
        #                         7GHZ      11GHZ      19GHZ     24GHZ     37GHZ     85GHZ   85GHZ     
        
    if band in [2,3,4]:
        coeff_a =  ant_approx_a[band]
        coeff_b =  ant_approx_b[band]
        coeff_c =  ant_approx_c[band]
        coeff_d =  ant_approx_d[band]
    else:
        raise ValueError(f'Band {band} not support for SSMIS')

    gain = coeff_a + coeff_b*np.exp(-coeff_c*delta) + np.exp(-coeff_d*delta*delta)
    return gain


def target_gain(delta_km,diameter_in_km = 30.0):
    import numpy as np

    #using the shape from AMSR2 so all shapes match

    coeff_a =  2.096e-6
    coeff_b =  4.059e-4
    coeff_c =  0.662
    coeff_d =  1.345

    scale_factor = 0.047933*30.0/diameter_in_km
    delta = delta_km*scale_factor
    gain = coeff_a + coeff_b*np.exp(-coeff_c*delta) + 	np.exp(-coeff_d*delta*delta)
    
    return gain
