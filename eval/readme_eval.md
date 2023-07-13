# This is description of the evaluation code in this folder (/eval)

## 1. Introduction

The purpose is to provide and estimate (or at least relative estimate) of the effect of footprnt shape and location errors on brightness temperatures for real-world scenes.  We use resampling to convert satellite measurements from a set of elliptical footprints at arbitrary locations to circular gaussian footprints on an Earth-fixed regular grid.  This is done in two steps -- resample to circular footprints in the satellite-swath geometry, and then either 1) Choose the closest such footprint to the desired location or 2) Interpolate the 4 surrounding footprints to the deired location.  Method 1 usually leads to some amount of location errors.  Method 2 results in an effective footprint that is a slightly stretched version of the circular target footprint.  An important question is which method results in the smallest radiance error for real world scenes.

For microwaves, some of the most challenging scenes are mixed land/water scenes.  The lower emissivity for water surfaces leads to a strong radiance difference between land and water.  This can range between 60K and 250K, depending on the frequency, polarization, surface temperature.  In the following analysis, we assume that this difference is 100K, a representative (and easy to remember) value.  This could be scaled to find a more accurate estimate by including information about the actual land/water contrast for a given channel and location.

## Approach

The approach is to use real-world land/water masks.  We focus on three regions, presented in order of how much challenge the spatial pattern presents to the resampling algorithm.  

1. Farmland in the upper midwest of the US (Iowa).  These scenes are mostly land, with isolated rivers and reservoirs.
2. A region south east of Hudson Bay in the Canadian Shield.  The geology in this area results in a large number of lakes.
3. A complex coastline from southern Maine.  There are numerous inlets and islands in this area.

The land water mask is dervied from the high resolution masks from Hansen et al, downsampled to 1 km resolution. 

There are three comparisons done.

1. "Ideal": Target footprint vs resampled footprint, both centered on the resampled location.
2. "Closest": Target footprint, translated to a random location in the quadrilateral surround the resampled footprint.  The quadrilateral is defined by 1/2 the distance to neighboring resampled locations, i.e. the locations where the resampled locations would be found to be the closest resampling location.  This is     what is currently done for the MWI Earth-grid resmapling.  In this case, there are mis-location errors of a up to 2-3 km.
3. "Interpolated":  The four closest resampled locations are interpolated to the translated target location.  Note that the shape of the resampled footprint is now distorted and no longer a nearly exact match to the traget footprint.

Translation of the target footprint involves re-locating it by arbitrary amounts in both the X and Y directions.  This is done using a spline-based interpolation routine from scipy. 

