#     From O:\ssmis\resampling\find_gains.f

#     These are Frank's comments on the fortran code.
#     !see memo2 in O:\ssmis\resampling

#  Comments below by CM, Summer 2022
#
#  converted to python so I can see what's going on
#  the python version is used to create the input data for 
#  the resampling weight calculations
#
#   real(8), parameter :: beamwidth(4)=(/1.87, 1.62, 1.05, 0.43/)  !see  'memo1.txt'
    #   subroutine fd_gain_targets(beamwidth,tht, gain)
    #   implicit none
 
 	#   real(8), parameter :: pi=3.141592653589793d0
	#   real(8), parameter :: rad=pi/180.d0
    #   real(8), parameter :: ln2=dlog(2.d0)
    #   real(8), parameter :: c1=4.d-4 
    #   real(8), parameter :: c2=0.456008940886013d0  !=sqrt(0.3d0*dlog(2.d0))

    #   real(8) beamwidth,tht, gain
	#   real(8) ththalf,ainv,x

	#   ththalf=0.5d0*beamwidth
    #   ainv=(rad*ththalf)**2/ln2
	#   x=tht/ththalf
	#   gain = (exp(-ln2*x*x) +  c1*exp(-c2*x))/(1.0021d0*pi*ainv) !1.0021 is fudge factor to bring norm to a little less than unity
    #   return
    #   end

 
def SSMI_antenna_gain(delta,freq):
    import numpy as np

    theta_full   = np.array([1.87, 1.62, 1.05, 0.43],dtype=np.float64)
    theta_half   = theta_full/2.0
    sin_anglim   = np.sin(3.0*theta_half*np.pi/180.0)
    ant_approx_a = np.array([2.00e-6,  2.00e-6,  2.00e-6,  2.00e-6,  2.00e-6],dtype=np.float64)
    ant_approx_b = np.array([4.00e-4,  4.00e-4,  4.00e-4,  4.00e-4,  4.00e-4],dtype=np.float64)
    
    ant_approx_c = np.zeros(4,dtype=np.float64)
    ant_approx_d = np.zeros(4,dtype=np.float64)

    ant_approx_d = -np.log(0.5)/theta_half**2
    
    ant_approx_c = np.sqrt(0.3*ant_approx_d)
    
    #   0         1         2         3 
    #   19GHZ     24GHZ     37GHZ     85GHZ     
        
    if freq in [0,1,2,3]:
        coeff_a =  ant_approx_a[freq]
        coeff_b =  ant_approx_b[freq]
        coeff_c =  ant_approx_c[freq]
        coeff_d =  ant_approx_d[freq]
    else:
        raise ValueError(f'Freq Index {freq} not valid for SSMI')

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

# Contents of memo1.txt from O:\ssmi\resampling - 3/15/2023
'''
        *** FINDING RESAMPLING WEIGHTS FOR SSMI ***
                                      
Here I find resampling weights for the SSMI.  There are two sets:
1.  The standard set I have been using that map 22 and 37 onto the 19 footprint
2.  A new set that regrids the obs onto both the A and B scans, all 128 samples.
Here I described the procedure and 'resampling.pptx' shows the resulting noise reduction factors and fit quality.

All looks good. 
 
My only concern is the relatively poor fit quality for regridding the 37 GHz onto even-numbered cells.
The Nyquist sampling along the scan (as compared to along the track) is not very good.
There is not way to avoid this, and I dont think it will cause any real-world problems.

As discussed below, the standard set of wts that map 22 and 37 onto the 19 footprint produce nearly identical results as the previous 
version.

The routine for apply the new wts to the TAs is 'O:\ssmi\routines\resampling_routines_l1c.f'
It has been thoroughly tested and is ready for ops processing.

During this analysis I had called it 'O:\ssmi\routines\resampling_routines_v8.f'
but on 3/15/2023 changed the name to 'O:\ssmi\routines\resampling_routines_l1c.f' 
I am using the suffix l1c to denote the new routines associated with the ssmi reprocessing effort.




==========================================================================================================
                                      *** F08 VERSUS THE OTHER SSMI ***
The analyses for the 5 ssmi other than f08 are done in this folder:  O:\ssmi\resampling
The analyses for  F08 are done in O:\ssmi\resampling_f8
The analyses are completely analogous.
I had to do a separate analysis for F08 because it looks backwards.
This matters at the scan  edges, where the neigboring scan order is reversed

I only had to change my reference orbit file.  Here I used F13 orbit 5000
For F08  is used F08 orbit 6000.
Some of the programs required no changes and the rest required modifying one or new lines.
The nobs_max was slightly different for the two analyses.
==========================================================================================================
                                      *** PREVIOUS ANALYSIS ***
                                      
The previous  derivation of resampling wts was done Aug 2006 in folder O:\ssmi5\V06\resampling,
It derived wts for mapping 18, 22, and 37 onto 18.  However the 18 onto 18 wts were never used.
The 2006 derviation used a slightly different model for the antenna, but the difs are too small to matter.
See 'new theoretical gain patterns' section below,
I verified the wts coming from this analysis produce essentially the same TAs.
See 'changes to standard 'o:\ssmi\routines\resample_ta_v8.f'  section below.
So this new analysis and wts completely  replace the 2006 results                                     
==========================================================================================================
                                      *** ANTENNA BEAMWIDTH  ***

I found a reference that gives the ssmi antenna beamwidth.  These are probably from Hollingers Report.
The reference is shown on the first page of'resampling.pptx'
The beamwidths are :
beamwidth(4)=(/1.87, 1.62, 1.05, 0.43/) 
The first 3 are the same as used in the 2006 derivation, which did not consider 85 ghz.
==========================================================================================================
                                 *** RUNNING THE PROGRAMS ***

The sequence of running the programs is the usual that I have done many times.
The programs and data files are in chronological order as they were run.
Look at he explorer listing for this folder to see the list of programs.
An analogous set of programs is on O:\ssmi\resampling_f8 for the f08 SSMI

==========================================================================================================

                          *** CHANGES TO STANDARD 'O:\SSMI\ROUTINES\RESAMPLE_TA_V8.F' ***
                          
1.  Shift to lo grid is now done in routine, since the program does not call convert_to_lores

2.  I found that the weight are only non-zero for +-4 scans and +-4 for cells about the cell be process.
    This is a much small range than the array allows for, which is +-18 and +-14 [resample_wt(ncel_rsp,-18:18,-14:14,2,2)]
    However, this is not really slowing the routine down by very much, so i did not make any mods.
    It does mean you get full resampling starting with scan 5 and ending with scan numscan-4

3.  Subroutine never use begin or end scans that are within 4 scans for end
    by doing this, jscan is never less that 1 or greater than numscan
    to make the change the following was done 
    if(iscan.lt.5 .or. iscan.gt.numscan_logrd-3) then  
    changed to 
    if(iscan.lt.5 .or. iscan.gt.numscan_logrd-4) then 
    
This new routine is called 'o:\ssmi\resampling\resample_ta_to_18ghz.f'
This new routine is tested by running 'test_resampling.f'.
This test program compares the ta coming from 'o:\ssmi\resampling\resample_ta_to_18ghz.f'
with the ta coming from  my new hi-grid routine 'o:\ssmi\routines\resampling_routines_v8.f'
The following listing shows the two sets of ta are virtually identical: rms dif < 0.01K.

O:\ssmi\resampling>test_resampling
 nbad_orbits =            1
 13  5000 1996  3 12  14.84
  0.000016  0.003142  0.000000  0.000000  0.000001  0.005240  0.000003  0.009069 -1.000000 -1.000000  0.098862  0.000000  0.094009  0.160706
 13  5001 1996  3 12  16.54
  0.000002  0.003278  0.000000  0.000000  0.000004  0.005386  0.000012  0.009267 -1.000000 -1.000000  0.074463  0.000000  0.082779  0.159958
norm end

This is an important result.  Is show:
1.  previous  resampling weights and the new ones derived herein produced the same tas.
    this is further confirmed by  'resampling.pptx'

2.  it verified the application procedure for apply the weight to the ta in ops processing

====================================================================================================================
                               *** NEW THEORETICAL GAIN PATTERNS ***
'ts_fd_gain_targets.f' shows that the theoretical gain given by the routine
'o:\mwi\antenna\fd_gain_targets.f' is nearly the same as that used in 
'O:\ssmi5\V06\resampling\find_gains.f' except that 
the 'ts_fd_gain_targets.f' has been normalized to be close to unity.
The difference in the normalization should have no effect on resampling weight.
When the off-boresight angle is greater than 1.5 full beamwidths, two expression diverage a little.
Run 'ts_fd_gain_targets.f' to see this.
Conclusion:  The two ways off computting the gain (old ssmi and mwi) will produce essentially the same resampling results.
             I will use the mwi routine.
=======================================================================================================================
                               *** ALONG-SCAN  SMEARING ***

	write(*,*) abs(dalpha)*scan_period/360. gives 8.440000000000001E-003

time between samples (19-37 ghz) = 8.440000000000001E-003
Reference to Hollinger's paper say integration times are 7.95 ms (19-37) and 3.89 ms (85)
So intergation time is 7.95/8.44 = 94% of the time
so there is 3% of dead time on either end
thus delta_xcel is -0.47 to + 0.47

	do  90 ismear=-5,5
	xcel=icel + 0.094*ismear

===========================================================================================
Setting smoothing factor from smooth_fac=1.d-18  to smooth_fac=0 have identical check_fit results.
This was suprising becasue i thought using  zero would produce a matrix inversion problem

===========================================================================================
                *** NON-ZERO RANGE OF RESAMPLING WEIGHT ***

The resampling wts have a range of -18:18. -14:14 (cells, scans).  This is typical for other sensors.
However, I found the wts were non-zero only over the range -4:4. -4:4.
I found this to be true for this new analysis as well.
I had not realized this before.  So in ops processing I only average over  +- 4 scans, not +-14.	

===========================================================================================
                  *** BEGIN AND END OF OORBITS
There are about 160 scan before and 160 afater an orbit, i.e. the ears.
Thus I can removed the first 10  scans (5 a/b pairs) and the last 8 scans (4 a/b pairs)
for which incomplete resampling occurs

'''
