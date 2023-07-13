import numpy as np
import xarray as xr
from rss_gridding.local_earth_grid import localearthgrid
#from hansen_land_fraction import get_hansen_land_fraction_local
#from common.read_gain_array import read_gain_array

import matplotlib.pyplot as plt
import pandas as pd
import os
print(os.getcwd())
print(__file__)
from common.random_land_mask import SomewhatRandomLandWater

class MakeRandomScene:

    def __init__(self,
                 scene='lakes',
                 subgrid_size=100,
                 beamwidth=70,
                 gradient_magnitude=0.5,
                 land_water_contrast=100.0,
                 scene_gen=None,
                 use_precomputed_locations=False,
                 num_random_samples=1000):
        
        self.scene = scene
        self.subgrid_size = subgrid_size
        self.beamwidth = beamwidth
        self.land_water_contrast = land_water_contrast
        self.gradient_magnitude = gradient_magnitude
        self.scene_gen = scene_gen
        self.use_precomputed_locations = use_precomputed_locations

        rng = np.random.default_rng()
        self._random_array = rng.random((num_random_samples,10))

        if scene in ['lakes','midwest','coastline']:
            self.scene_gen = SomewhatRandomLandWater(scene=scene)

        if self.use_precomputed_locations:
            path_to_random_locs = 'L:/access/resampling/SSMI/f13/eval/random_locations'
            txt_file = path_to_random_locs / f'random_locations_{scene}_{60}km.txt'
            df = pd.read_csv(txt_file,header=None)
            lons = df.iloc[:,0].values
            lats = df.iloc[:,1].values
            xshft = df.iloc[:,2].values
            yshft = df.iloc[:,3].values
            if scene in ['edges','gradient']:
                phis = df.iloc[:,4].values

            self.random_lons = lons
            self.random_lats = lats
            self.random_xshft = xshft
            self.random_yshft = yshft
            if scene in ['edges','gradient']:
                self.phis = phis

    def make_random_scene(self,*,
                        isamp):

        if self.scene in ['lakes','midwest','coastline']:
            grid,lat,lon = self.scene_gen.get_somewhat_random_earth_scene(subgrid_size=self.subgrid_size)

        elif self.scene == 'edges':
            # make grid with edge

            # for the edge case, we want to make a grid with a single edge
            # it needs to be at a random angle, and a random distance 
            # from the center.  The distance is chosen to be within 1/2 of
            # the beamwidth.  The angle is chosen randomly.

            lon = 180.0
            lat = 0.0
            grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                            delta_x = 1.,delta_y = 1.,
                            num_x = 1+2*self.subgrid_size,num_y = 1+2*self.subgrid_size,
                            name='Land Fraction')
            
            # chose a random angle for the edge
            if self.use_precomputed_locations:
                phi = self.phis[isamp]
                offsetx = self.random_xshft[isamp]
                offsety = self.random_yshft[isamp]
            else:
                phi = 2.0*np.pi*self._random_array[isamp,0]
                offsetx = -0.5*self.beamwidth + self.beamwidth*self._random_array[isamp,1]
                offsety = -0.5*self.beamwidth + self.beamwidth*self._random_array[isamp,2]
            
            
            #phi_all.append(phi)
            dx = np.cos(phi)
            dy = np.sin(phi)

            edge = np.zeros((2*self.subgrid_size+1,2*self.subgrid_size+1))
            for ix in range(-self.subgrid_size,self.subgrid_size+1):
                for iy in range(-self.subgrid_size,self.subgrid_size+1):
                    if ((dx*(ix-offsetx))+(dy*(iy-offsety))) > 0.0:
                        edge[ix+self.subgrid_size,iy+self.subgrid_size] = 1.0
            grid.setdata(edge,'Land Fraction')

        elif self.scene == 'gradient':
            # make grid spatial gradient

            # for the gradient case, the gradient is along a random angle
            # translation is not required

            lon = 180.0
            lat = 0.0
            grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                            delta_x = 1.,delta_y = 1.,
                            num_x = 1+2*subgrid_size,num_y = 1+2*subgrid_size,
                            name='Land Fraction')

            if self.use_precomputed_locations:
                phi = self.phis[isamp]
            else:
                phi = 2.0*np.pi*self._random_array[isamp,0]

            #phi_all.append(phi)
            dx = np.cos(phi)
            dy = np.sin(phi)
            grad = np.zeros((2*subgrid_size+1,2*subgrid_size+1))
            for ix in range(-subgrid_size,subgrid_size+1):
                for iy in range(-subgrid_size,subgrid_size+1):
                    dot = dx*ix +dy*iy
                    grad[ix+subgrid_size,iy+subgrid_size] = (self.gradient_magnitude/self.land_water_contrast)*dot
            grid.setdata(grad,'Land Fraction')
        else:
            raise ValueError(f'Scene: {scene} not implemented')
    
        return grid,lat,lon

if __name__ == "__main__":

    print('Testing MakeRandomScene')
    for scene in ['lakes','midwest','coastline','edges','gradient']:

        scene_gen = MakeRandomScene(scene=scene,subgrid_size=100,beamwidth=70)

        isamp = 1
        grid,lat,lon = scene_gen.make_random_scene(scene=scene,subgrid_size=100,beamwidth=70,isamp=isamp)

        grid.contourplot(vmin=-1.0,vmax=1.0,plt_contours=False)

    print