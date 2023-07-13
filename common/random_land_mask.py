import sys

sys.path.append("M:/job_access/python/land_water/")
sys.path.append("M:/job_access/python/resample_wts/")

import numpy as np
from numpy.random import default_rng
from scipy.interpolate import RectBivariateSpline

#local
from hansen_land_fraction import get_hansen_land_fraction_local,read_4_hansen_land_fraction_file_that_contains
from rss_plotting.plot_2d_array import plot_2d_array
from rss_gridding.local_earth_grid import localearthgrid

class SomewhatRandomLandWater:

    def __init__(self,
                 scene : str=None,
                 center_lat=None,
                 center_lon=None,
                 num_samples_possible = 10000):

        #get a large random array        
        rng = default_rng()
        self.random_array = rng.random((num_samples_possible,10))
        self.sample_index = 0
        self.num_samples_possible = num_samples_possible

        #set up local area map
        if scene is not None:
            if scene in ['lakes','midwest','coastline']:
                if scene == 'lakes':
                    center_lat = 53.0
                    center_lon = 295.0
                    dlatlon = 2.0
                elif scene == 'midwest':
                    center_lat = 45.2
                    center_lon = 262.0
                    dlatlon = 2.0
                elif scene == 'coastline':
                    center_lat = 43.5
                    center_lon = 290.  
                    dlatlon = 1.0
                else:
                    raise ValueError(f'Scene: {scene} not implemented')
            
        if ((center_lat is  None) or (center_lon is None)):
            raise ValueError('center_lat and center_lon must be provides if scene is not')

        self.center_lon = center_lon
        self.center_lat = center_lat
        self.dlatlon = dlatlon

        lf_dataset = read_4_hansen_land_fraction_file_that_contains(lon=center_lon,lat=center_lat,verbose = True)
        
        self.land_fraction = lf_dataset['land_fraction']
        self.latitudes = lf_dataset['Latitude']
        self.longitudes = lf_dataset['Longitude']

        self.interp_spline = RectBivariateSpline(self.latitudes,self.longitudes,self.land_fraction)

        fig,ax = plot_2d_array(self.land_fraction,self.longitudes,self.latitudes,zrange=[0.0,1.0],cmap='Blues_r')
        lons_to_plot = np.array([self.center_lon-self.dlatlon,
                        self.center_lon+self.dlatlon,
                        self.center_lon+self.dlatlon,
                        self.center_lon-self.dlatlon,
                        self.center_lon-self.dlatlon])
        lats_to_plot = np.array([self.center_lat-self.dlatlon,
                        self.center_lat-self.dlatlon,
                        self.center_lat+self.dlatlon,
                        self.center_lat+self.dlatlon,
                        self.center_lat-self.dlatlon])    
        
        lons_to_plot = (lons_to_plot+360.0+180.0)%360.0 - 180
        ax.plot(lons_to_plot,lats_to_plot,color='red')
        self.fig = fig
        self.ax = ax

    def get_local_lf(self,lon=262.23,lat=52.43,size_km = 50,verbose = True):

        grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 1+2*size_km,num_y = 1+2*size_km,
                          name='Land Fraction')

        lats = grid.lats
        lons = grid.lons
        lf_on_local_grid = self.interp_spline(lats,lons,grid=False)

        grid.setdata(lf_on_local_grid,name='Land Fraction')

        return grid

    def get_somewhat_random_earth_scene(self,
                                        subgrid_size=70,
                                        lat_in = None,
                                        lon_in = None):

        # if lat_in and lon_in are provided, use them
        # otherwise, use the random array to get a new lat/lon within
        # 0.5 degree of the center lat/lon

        if lat_in is None:
            dlat = self.dlatlon*(-0.5 + 1.0*self.random_array[self.sample_index,1])
            lat = self.center_lat + dlat
        else:
            lat = lat_in

        if lon_in is None:
            dlon = self.dlatlon*(-0.5 + 1.0*self.random_array[self.sample_index,0])
            lon = self.center_lon + dlon
        else:
            lon = lon_in

        self.sample_index = (self.sample_index + 1) % self.num_samples_possible

        grid = self.get_local_lf(lon=lon,lat=lat,size_km = subgrid_size,verbose=False)
        return grid,lat,lon

if __name__ == '__main__':

    import matplotlib.pyplot as plt

    lakes_scene_gen = SomewhatRandomLandWater(scene='lakes')

    lf_grid = lakes_scene_gen.get_somewhat_random_earth_scene()

    fig,ax,im = lf_grid.contourplot(name = 'Land Fraction',plt_contours=False)

    lf_grid = lakes_scene_gen.get_somewhat_random_earth_scene()

    fig,ax,im = lf_grid.contourplot(name = 'Land Fraction',plt_contours=False)

    lf_grid = lakes_scene_gen.get_somewhat_random_earth_scene(subgrid_size=100)

    fig,ax,im = lf_grid.contourplot(name = 'Land Fraction',plt_contours=False)

    coast_scene_gen = SomewhatRandomLandWater(scene='coastline')
    lf_grid = coast_scene_gen.get_somewhat_random_earth_scene(subgrid_size=100)
    fig,ax,im = lf_grid.contourplot(name = 'Land Fraction',plt_contours=False)
    plt.show()
    print