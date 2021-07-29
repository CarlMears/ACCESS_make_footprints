Each satellite subdirectory contains the following routines, at a minimum

*sat*_Antenna_Gain.py
    This routine provides the antenna gain as a function of off-boresight angle.

    *sat*_Antenna_Gain(delta,band)  delta is the off-boresight angle in degrees
                                    band is the band number, which is instrument specific to some degree, but
                                    in general, the assignments should look like this..
                                    0 is 6/7 GHz
                                    1 is 11 GHz
                                    2 is 19 GHz
                                    3 is 23/24 GHz
                                    4 is 36/37 GHz
                                    5 is 85/89/91 GHz

make_*sat*_footprints_for_resampling.py
    This routine reads some sort of file with footprint locations, and then calculates target and source
    patterns at the appropriate subset of those locations.
    

Once these are set up, the actual calculation of the fooprint weights is done in fortran
