import numpy as np
from common.LocalEarthGrid import LocalEarthGrid
import xarray as xr
from rss_gridding.local_earth_grid import localearthgrid
from hansen_land_fraction import get_hansen_land_fraction_local
from common.read_gain_array import read_gain_array

import matplotlib.pyplot as plt
from numpy.random import RandomState
from numpy.random import default_rng

from common.random_land_mask import SomewhatRandomLandWater


def make_random_scene(*,scene,scene_gen=None,subgrid_size):

    if scene in ['lakes','midwest','coastline']:
        grid,lat,lon = scene_gen.get_somewhat_random_earth_scene(subgrid_size=subgrid_size)

    elif scene == 'edges':
        # make grid with edge

        # for the edge case, we want to make a grid with a single edge
        # it needs to be at a random angle, and a random distance 
        # from the center.  The distance is chosen to be within 1/2 of
        # the beamwidth.  The angle is chosen randomly.

        lon = 180.0
        lat = 0.0
        grid = LocalEarthGrid(center_lon = lon,center_lat = lat,
                        delta_x = 1.,delta_y = 1.,
                        num_x = 1+2*subgrid_size,num_y = 1+2*subgrid_size,
                        name='Land Fraction')
        
        # chose a random angle for the edge
        if use_precomputed_locations:
                phi = phis[isamp]
        else:
            phi = 2.0*np.pi*random_array[isamp,0]

        # chose a random offset from the center
        offsetx = -0.5*beamwidth + beamwidth*random_array[isamp,1]
        offsety = -0.5*beamwidth + beamwidth*random_array[isamp,2]
        
        
        phi_all.append(phi)
        dx = np.cos(phi)
        dy = np.sin(phi)

        edge = np.zeros((2*subgrid_size+1,2*subgrid_size+1))
        for ix in range(-subgrid_size,subgrid_size+1):
            for iy in range(-subgrid_size,subgrid_size+1):
                if ((dx*(ix-offsetx))+(dy*(iy-offsety))) > 0.0:
                    edge[ix+subgrid_size,iy+subgrid_size] = 1.0
        grid.setdata(edge,'Land Fraction')

    elif scene == 'gradient':
        # make grid spatial gradient

        # for the gradient case, the gradient is along a random angle
        # translation is not required

        lon = 180.0
        lat = 0.0
        grid = LocalEarthGrid(center_lon = lon,center_lat = lat,
                        delta_x = 1.,delta_y = 1.,
                        num_x = 1+2*subgrid_size,num_y = 1+2*subgrid_size,
                        name='Land Fraction')

        if use_precomputed_locations:
            phi = phis[isamp]
        else:
            phi = 2.0*np.pi*random_array[isamp,0]

        phi_all.append(phi)
        dx = np.cos(phi)
        dy = np.sin(phi)
        grad = np.zeros((2*subgrid_size+1,2*subgrid_size+1))
        for ix in range(-subgrid_size,subgrid_size+1):
            for iy in range(-subgrid_size,subgrid_size+1):
                dot = dx*ix +dy*iy
                grad[ix+subgrid_size,iy+subgrid_size] = (gradient_magnitude/land_water_contrast)*dot
        grid.setdata(grad,'Land Fraction')
    else:
        raise ValueError(f'Scene: {scene} not implemented')
