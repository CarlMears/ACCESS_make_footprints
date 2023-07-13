


import sys

sys.path.append("M:/job_access/python/land_water/")
sys.path.append("M:/job_access/python/resample_wts/")
import numpy as np
import xarray as xr
from rss_gridding.local_earth_grid import localearthgrid
from hansen_land_fraction import get_hansen_land_fraction_local
from common.read_gain_array import read_gain_array

import matplotlib.pyplot as plt
from numpy.random import RandomState
from numpy.random import default_rng

def circle(*,diameter):
    phi = np.arange(0,2.0*np.pi,0.001)
    x = diameter/2.0 * np.cos(phi)
    y = diameter/2.0 * np.sin(phi)
    return (x,y)

rng = default_rng()
scene_list = ['lakes','midwest','coastline','linear','gradient']
panel_labels = ['a','b','c','d','e']

fig,axs = plt.subplots(ncols = 5,nrows=1,figsize=(20,4))
for scene_index,scene in enumerate(scene_list):
    l = 50
    if scene in ['lakes','midwest','coastline']:
        if scene == 'lakes':
            center_lat = 52.0
            center_lon = 296.0
        elif scene == 'midwest':
            center_lat = 42.9
            center_lon = 262.0
        elif scene == 'coastline':
            center_lat = 43.7
            center_lon = 295.0
        else:
            raise ValueError(f'Scene: {scene} not implemented')

        grid = get_hansen_land_fraction_local(lon=center_lon,lat=center_lat,size_km = 50,verbose=False)
    elif scene == 'linear':
        # make grid with edge
        size_km = 50
        lon = 180.0
        lat = 0.0
        grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                        delta_x = 1.,delta_y = 1.,
                        num_x = 1+2*size_km,num_y = 1+2*size_km,
                        name='Land Fraction')
        phi = 2.0*np.pi*0.243
        
        offsetx = -5.0 + 10.0*0.421
        offsety = -5.0 + 10.0*0.789
        dx = np.cos(phi)
        dy = np.sin(phi)
        edge = np.zeros((101,101))
        for ix in range(-50,51):
            for iy in range(-50,51):
                if ((dx*(ix-offsetx))+(dy*(iy-offsety))) > 0.0:
                    edge[ix+50,iy+50] = 1.0
        grid.setdata(edge,'Land Fraction')
        grid.addlayer(name='Tb')
        tb_edge = 150.0+100.0*edge
        grid.setdata(tb_edge,'Tb')
    elif scene == 'gradient':
        # make grid with edge
        size_km = 50
        lon = 180.0
        lat = 0.0
        grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                        delta_x = 1.,delta_y = 1.,
                        num_x = 1+2*size_km,num_y = 1+2*size_km,
                        name='Land Fraction')
        phi = 2.0*np.pi*0.569
        
        dx = np.cos(phi)
        dy = np.sin(phi)
        edge = np.zeros((101,101))
        for ix in range(-50,51):
            for iy in range(-50,51):
                dot = dx*ix +dy*iy
                edge[ix+50,iy+50] = 0.5 + 0.01*dot
        grid.setdata(edge,'Land Fraction')
        grid.addlayer(name='Tb')
        grid.setdata(150.0+100.0*edge,'Tb')
    else:
        raise ValueError(f'Scene: {scene} not implemented')
    
    plt.rcParams.update({'font.size': 18})
    from matplotlib import cm
    cmap = cm.winter
    
    fig,axs[scene_index],im = grid.contourplot(name='Tb',vmin=150.0,vmax =250.0,fig_in=fig,
                         ax_in=axs[scene_index],
                          cbticks = [150.0,175.0,200.0,225.0,250.0],cmap = cmap,scale=False,
                          units = 'Brightness Temperature',plt_contours=False,plt_colorbar=False,
                          xlabel='EW distance (km)',ylabel='NS distance (km)')
    circ=circle(diameter = 30)
    axs[scene_index].plot(circ[0],circ[1],color='red')
    axs[scene_index].set_title(scene)
    axs[scene_index].text(-46,40,panel_labels[scene_index],color='red')

fig.subplots_adjust(wspace=0.3)

cb_ax = fig.add_axes([0.92,0.18,0.01,0.62])
cbar = fig.colorbar(im,cax=cb_ax,orientation='vertical')
cbar.set_label('Brightness Temperature (K)',fontsize=12)
cb_ax.tick_params(axis='y', labelsize=11)   

png_file = 'M:/job_access/docs/ResamplingPaper/figures/Figure_06/Figure_06.png'
tif_file = 'M:/job_access/docs/ResamplingPaper/figures/Figure_06/Figure_06.tif'
fig.savefig(png_file,bbox_inches='tight')
fig.savefig(tif_file,bbox_inches='tight')


                
plt.show()
       
print()
