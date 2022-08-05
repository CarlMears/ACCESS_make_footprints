Each satellite subdirectory contains the following routines, at a minimum

*sat*_Antenna_Gain.py

This routine provides the antenna gain as a function of off-boresight angle.

## *sat*_Antenna_Gain(delta,band)  
* delta is the off-boresight angle in degrees
* band is the band number, which is instrument specific to some degree, but
in general, the assignments should look like this..
- 0 is 6/7 GHz
- 1 is 11 GHz
- 2 is 19 GHz
- 3 is 23/24 GHz
- 4 is 36/37 GHz
- 5 is 85/89/91 GHz

## make_*sat*_footprints_for_resampling.py
This routine reads some sort of file with footprint locations, and then calculates target and source patterns at the appropriate subset of those locations.

The reason that I do it this way is that I don't need all the geolocation machinery -- I can just
use the locations that were already calculated for
the L1 data for the various satellites.  I just need to figure out where the swath crosses the equator
    
Once these are set up, the actual calculation of the fooprint weights is done in fortran.  This is set up to run in linux on the VM called CentOS-MWI-sim-env on CMEARS.  The active code in is in /opt/ACCESS/resampling/src/earth_grid

The Linux calculation generates two types of output.  First, it generates a big file of resample weights - e.g.

L:\access\resampling\AMSR2\resample_weights\circular_30\resample_weights_band_02.dat

It also generates files with resmapled gain patterns for each resmapled locations - e.g. for scan 01, cell 11:

L:\access\resampling\AMSR2\resample_gains_v3\circular_30km\band_02\s01c0011.dat

These second files can be plotted and used for other analysis.








