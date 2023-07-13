
import numpy as np
import matplotlib.pyplot as plt
import struct
import os
import sys
sys.path.append('M:/job_access/python/resample_wts/common/')
from read_gain_array import read_gain_array

def plot_2d_array(a, xvals, yvals,zrange=(0.0,1.0), title='', xtitle='', ytitle='',cmap='BrBG',plt_colorbar = True,zlabel=' ',fig=None,subplot=111):

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import copy

    X, Y = np.meshgrid(xvals, yvals)
    if fig is None:
        fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(subplot, title=title, xlabel=xtitle, ylabel=ytitle)
    im = ax.pcolormesh(X, Y, a, cmap=cmap, norm=colors.Normalize(vmin=zrange[0], vmax=zrange[1]))
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(20)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(16)

    if plt_colorbar:
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel(zlabel, fontsize=16)
    return fig,ax



if __name__ == "__main__":

    sat_name = 'SSMI'
    beamwidth = 70
    freq_num = 0

    gain_type='target'
    fovs_to_show =[1,16,32]
    kscan = 1
    if gain_type == 'source':
        fovs_to_show = list(np.arange(3,63,10))
        kscan=14
    total = np.zeros((1602,2402))
    for fov in fovs_to_show:
        target_gain,ilat_max,ilon_max = read_gain_array(sat_name=sat_name,
                                                        ksat=13,
                                                        beamwidth=beamwidth,
                                                        band_num=freq_num,
                                                        kscan=kscan,
                                                        fov=fov,
                                                        gain_type=gain_type,
                                                        verbose = True)

        total += target_gain
        print(fov,ilat_max,ilon_max)

    
    x_vals = np.arange(0,2402)
    y_vals = np.arange(0,1602)

    fig = plt.figure(figsize=(18,8))
    plot_2d_array(total,x_vals,y_vals,fig=fig)
    plt.show()

    print