
"""
   Calculations on a Local Earth Grid derived from a transverse Mercator projection
   
   Author:
       Carl A. Mears, Remote Sensing Systems
       
"""
import math
from math import sin,cos
import numpy as np
import collections
import pyproj  # does various projections
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import struct

class LocalEarthGrid():
    '''
    Class for calculating and displaying data on local cartesion earth grid
    
    Projection is transverse Mercator, centered at the center of grid.
    When grid is set up, lat and lon of each grid point is calculated
    Any number of layers can be added that share the common grid scheme
    
    Args:
        center_lon:    longitude (degrees, 0.0..360.0) of center of local grid
        center_lat:    latitude (degrees, -90.0...90.0) of center of local grid
        delta_x:       horizonal (East-West) grid spacing in km
        delta_y:       vertical (North-South) grid spacing in km
        num_x:         integer number of grid locations in horizontal direction
        num_y:         integer number of grid locations in vertical direction
                       numx and numy are usually odd, so that a grid point exists at
                       the exact center of the grid
        name:          string containing the name of the first layer in the class instance
       
       '''
       
    def __init__(self,center_lon=0.0,center_lat = 0.0,
                 delta_x = 2.0,delta_y = 2.0,
                 num_x = 201,num_y = 201,name='Default',constant_value = 0.0):
        '''Initialize the grid instance centered at center_lon,center_lat.  
           delta_x and delta_y are the grid spacing in km.  num_x and num_y 
           are the number of grid points in the x and y directions.  name is
           the name of the first layer in the grid instance.  constant_value
           is the value to which the first layer is initialized.'''
        self.center_lon = center_lon
        self.center_lat = center_lat
        self.delta_x = delta_x       # grid spacing in x
        self.delta_y = delta_y       # grid spacing in y
        
        if num_x%2 == 0:
            num_x = num_x+1
        if num_y%2 == 0:
            num_y = num_y+1
        
        self.num_x = num_x        # number of x grid points
        self.num_y = num_y        # number of y grid points
        
        self.x_extent = [-1*((self.num_x-1)/2)*self.delta_x,(self.num_x-1)/2*self.delta_x]
        self.y_extent = [-1*((self.num_y-1)/2)*self.delta_y,(self.num_y-1)/2*self.delta_y]
    
        ## set up array with the data
        self.data = collections.OrderedDict()
        self.data[name] = constant_value*np.ones([self.num_y,self.num_x])   ## set up first map with zeros
        
        ## popluate arrays with locations in km
        x_vals = np.arange(self.x_extent[0],self.x_extent[1]+self.delta_x,self.delta_x)        
        y_vals = np.arange(self.y_extent[0],self.y_extent[1]+self.delta_y,self.delta_y)
        self.xlocs,self.ylocs = np.meshgrid(x_vals,y_vals)
        
        # compute lats and lons using transverse mercator projection centered at center_lon and canter_lat
        # first, set up projection
        
        proj_command = '+proj=tmerc +lat_0='+str(self.center_lat) +' +lon_0=' + str(self.center_lon)
        proj_command = proj_command + ' +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
        transverse_mercator = pyproj.Proj(proj_command) 
        self.projection = transverse_mercator
        # do the transformation
        self.lons,self.lats = transverse_mercator(1000.0*self.xlocs,1000.0*self.ylocs,inverse=True)
        #self.lons and self.lats are arrays of locations in lat and lon corresponding to the 
        #locations in x,y space.
        
    def addlayer(self,name='Default'):
        '''Add another layer to the grid instance.'''
        if name in self.data:
            pass
            #print('Layer with name ='+name+' already exists - no layer added')
        else:
            self.data[name] = np.zeros([self.num_y,self.num_x])
            
    def haslayer(self,name='Default'):
        '''Check for the existance of a layer called name'''
        if name in self.data:
            return True
        else:
            return False
            
    def gridinfo(self):
        '''Return human-readable info about the grid parameters'''
        
        s =     'Center Longitude:  '+str(self.center_lon) + '\n'
        s = s + 'Center Latitude:     '+str(self.center_lat) + '\n'
        s = s + 'Number X Dimension: '+str(self.num_x) + '\n'
        s = s + 'Number Y Dimension: '+str(self.num_y) + '\n'
        s = s + 'Delta X:            '+str(self.delta_x) + '\n'
        s = s + 'Delta Y:            '+str(self.delta_y) + '\n'
        return s
        
    def addgaussian_at_lonlat(self,lon,lat,major,minor,azimuth,amplitude,normalize = False,name = None):
        '''
        Add a gaussian footprint to the specified layer
        
        Args:
            lon:       longitude of footprint center 
            lat:       latitude of footprint center 
            major:     major axis of 3db edge of footprint in km
            minor:     minor axis of 3 db edge of footprint in km
            azimith:   azimuth angle of footprint.  0.0 means major axis aligned East-West
            amplitude: overall amplitude of the footprint center
            normalize: if True, the footprint is normalized (numerically) so the the total over the
                       grid is equal to 1.0
            name:      the footprint is added to the layer with this name
            
        Returns:
            nothing
            
        '''
        
       
        if name is None:
            name = self.data.keys([0])
        ## find x and y for this longitude and latitude
        
        x,y = self.projection(lon,lat)
        x = x/1000.0
        y = y/1000.0
        
        sigma_x = major/2.355
        sigma_y = minor/2.355
        
        cos_azimuth = math.cos((math.pi/180.0)*azimuth)
        sin_azimuth = math.sin((math.pi/180.0)*azimuth)
        sin_2azimuth = math.sin(2.0*(math.pi/180.0)*azimuth)
    
        a = ((cos_azimuth*cos_azimuth)/(2.0*sigma_x*sigma_x) + 
             (sin_azimuth*sin_azimuth)/(2.0*sigma_y*sigma_y))
        b = (-1.0*sin_2azimuth/(4.0*sigma_x*sigma_x) + 
                  sin_2azimuth/(4.0*sigma_y*sigma_y))
        c = ((sin_azimuth*sin_azimuth)/(2.0*sigma_x*sigma_x) + 
             (cos_azimuth*cos_azimuth)/(2.0*sigma_y*sigma_y))

        gaussian = amplitude*np.exp(-1.0*(a*(self.xlocs-x)*(self.xlocs-x)+ 
                                     2.0*b*(self.xlocs-x)*(self.ylocs-y)+ 
                                         c*(self.ylocs-y)*(self.ylocs-y)))
        if normalize:
            tot = np.sum(gaussian)
            gaussian = gaussian/tot
        self.data[name] = np.add(self.data[name],gaussian)
    
    def getdata(self,name = None):
        '''Returns a deep copy of data from a layer as a numpy array.'''
        if name is None:
            name = self.data.keys([0])
            
        dataout = np.zeros([self.num_y,self.num_x])
        np.copyto(dataout,self.data[name])
        return dataout
        
    def setdata(self,newdata,name=None):
        '''Set the data in a layer to the supplied numpy array.'''
        if name is None:
            name = self.data.keys([0])
            
        if isinstance(newdata,np.ndarray):
            shp = np.shape(newdata)
            if ((shp[1] == self.num_x) and (shp[0] == self.num_y)):
                np.copyto(self.data[name],newdata)
            else:
                print('Array size does not match -- no assignment made')
        else:
            print('Argument is not a numpy array -- no assignment made')

    def adddata(self,newdata,name=None):
        '''Set the data in a layer to the supplied numpy array.'''
        if name is None:
            name = list(self.data.keys())[0]
            
        if isinstance(newdata,np.ndarray):
            shp = np.shape(newdata)
            if ((shp[1] == self.num_x) and (shp[0] == self.num_y)):
                np.copyto(self.data[name],newdata+self.data[name])
            else:
                print('Array size does not match -- no addition performed')
        else:
            print('Argument is not a numpy array -- no addition performed')
    
    def zero(self,name='None'):
        '''Set the data in a layer to zero.'''
        self.setdata(np.zeros((self.num_x,self.num_y)),name)
            
    def compare(x,y):
        '''Compare two grid instances to see if they are compatible'''
        if ((not isinstance(x,LocalEarthGrid)) or 
            (not isinstance(y,LocalEarthGrid))): 
                return False
        
        grids_same = True
        if x.__num_x != y.__num_x:
            grids_same = False
        if x.__num_y != y.__num_y:
            grids_same = False
        if x.__delta_x != y.__delta_x:
            grids_same = False         
        if x.__delta_y != y.__delta_y:
            grids_same = False   
        if x.__center_lat != y.__center_lat:
            grids_same = False
        if x.__center_lon != y.__center_lon:
            grids_same = False
            
        return grids_same
        
        
    def multiply(self,name1,name2):
        '''Multiply data from two layers and return product as numpy array.'''
        new_grid = np.copy(self.getdata(name1)) * np.copy(self.getdata(name2))
        return new_grid
        
    def add(self,name1,name2):
        '''Add data from two layers and return product as numpy array.'''
        new_grid = np.copy(self.getdata(name1)) + np.copy(self.getdata(name2))
        return new_grid

    def subtract(self,name1,name2):
        '''Difference data from two layers and return product as numpy array.'''
        new_grid = np.copy(self.getdata(name1)) - np.copy(self.getdata(name2))
        return new_grid
        
    def plot3d(self,name=None):
        '''Make a 3-D plot of the data in layer.'''

        
        if name is None:
            name = self.data.keys([0])
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        surf = ax.plot_surface(self.xlocs, self.ylocs, self.data[name], 
                               rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(-1.01, 3.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()

    def minmax(self,name):
        data = self.data[name]
        maxval = np.nanmax(data)
        minval = np.nanmin(data)
        return minval,maxval
        
    def total(self,name,no_area = False):
        '''Compute total of all data in layer'''
        data = self.data[name]
        area_per_pixel = self.delta_x *self.delta_y
        if no_area:
            area_per_pixel = 1.0
        total = area_per_pixel * np.sum(data)
        return total
    
    def mean_position(self,name):
        '''Compute the mean location of data in layer, if taken as weights'''
        x_position = np.sum(self.xlocs * self.data[name])/np.sum(self.data[name])
        y_position = np.sum(self.ylocs * self.data[name])/np.sum(self.data[name])
        return x_position,y_position
        
    def scale(self,name,scale_factor=1.0):
        '''Scale data in a layer by a scalar scale factor'''
        self.data[name] = np.copy(self.data[name])*scale_factor
        return 

    def normalize(self,name):
        tot = np.sum(self.data[name])
        self.data[name] = self.data[name]/tot
        
    def calc_difference_stats(self,name1,name2,normalize = False,verbose = False):
        data1 = self.data[name1]
        data2 = self.data[name2]
        area_per_pixel = self.delta_x *self.delta_y
        total1 = area_per_pixel * np.sum(data1)
        if verbose:
            print('Total Area '+name1+' '+str(total1))
        total2 = area_per_pixel * np.sum(data2)
        if verbose:
            print('Total Area '+name2+' '+str(total2))
        if normalize:
            data2 = data2*total1/total2
        diff = data1 - data2
        sum_sqr_diff = area_per_pixel*np.sum(diff * diff)
        if verbose:
            print('Sum Squared Diff '+str(sum_sqr_diff))
        return sum_sqr_diff
    
    def plot_ellipse(self,x0,y0,a,b,phi,color='blue',alpha = 1.0,facecolor = 'red',fig='None',sbplt = (1,1,1),ax=None):
        '''Plot an ellipse on the specified subplot.'''
        
        pi = 3.14159265359
        t = 0
        x = np.zeros([81])
        y = np.zeros([81])
        phi = 2.0*pi*phi/360.0
        for i in range(0,81):
            t = i*2.0*pi/40.0    
            x[i] = x0+a*cos(t)*sin(phi)-b*sin(t)*cos(phi)
            y[i] = y0+a*cos(t)*cos(phi)+b*sin(t)*sin(phi)
            
        if ax is None:
            ax = fig.add_subplot(sbplt[0],sbplt[1],sbplt[2])
        if facecolor is None:
            ax.plot(x,y,color=color,alpha=alpha)
        else:
            ax.fill(x,y,edgecolor=color,alpha=alpha,facecolor = facecolor)
        ax.set_xlim(self.x_extent[0],self.x_extent[1])
        ax.set_ylim(self.y_extent[0],self.y_extent[1])
        
        return ax

    def write_binary(self,outfile = None,name=None,save_thres=0.00001):



        if name is None:
            name = self.data.keys([0])

        z = self.data[name]
        good_loc = np.where(z > save_thres)

        with open(outfile,'wb') as f:
            for i in np.arange(0,good_loc[0].shape[0]):
                    iy = good_loc[0][i]
                    ix = good_loc[1][i]
                    z2 = z[iy,ix]
                    if z2 < 0.0000099:
                        raise ValueError()
                    s = struct.pack('<hhd',iy,ix,z2)
                    f.write(s)

    def print_amplitude_stats(self,name=None):
        if name is None:
            name = self.data.keys([0])

        z = self.data[name]
        print(f'Max value = {np.max(z)}')
        for scale_factor in [0.1,0.01,0.001,0.0001,0.00001]:
            print(f'Scale Fact {scale_factor}: {np.sum(z > scale_factor*np.max(z))}')


    def contourplot(self,outfile = None,name=None,vmin=-1.1,vmax =1.1,fig_in=None,
                    ax_in=None,sbplt=[1,1,1],cbticks = [-1.0,-0.5,0.0,0.5,1.0],cmap = None,
                    nodata=False,title = None,scale=True,scale_factor=None,units=' ',
                    plt_contours=True,plt_colorbar=True,
                    xlabel='EW distance (km)',ylabel='NS distance (km)',
                    panel_label=None):
        '''Make contour plot of data in specified layer, on the specified subplot.'''
        
        if name is None:
            name = list(self.data.keys())[0]
        if fig_in is None:
            fig = plt.figure(figsize = (9,9))
        else:
            fig = fig_in
        if cmap is None:
            cmap = cm.BrBG
            
        if ax_in is None:
            ax = fig.add_subplot(sbplt[0],sbplt[1],sbplt[2])
        else:
            ax = ax_in
        #plt.clf()
        ax.tick_params(axis='x')
        ax.tick_params(axis='y')
        #data = self.data[name]
        data = self.getdata(name = name)
        if scale is True:
            if scale_factor is None:
            # rescale the data to have 1.0 as maximum
                datamax = np.max(data)
                if datamax > 0:
                    data = data/np.max(data)
                else:
                    data = data/np.min(data)
            else:
                data = data/scale_factor

        if not nodata:
            im = ax.imshow(data, interpolation='bilinear', origin='lower',
                           cmap=cmap,extent=(self.x_extent[0],
                                             self.x_extent[1],
                                             self.y_extent[0],
                                             self.y_extent[1]),
                           norm =  plt.Normalize(vmin=vmin, vmax=vmax))
                    #vmin = vmin,vmax = vmax)
        else:
            im = ax.imshow(data, interpolation='bilinear', origin='lower',
                           cmap=cmap,extent=(self.x_extent[0],
                                             self.x_extent[1],
                                             self.y_extent[0],
                                             self.y_extent[1]),
                        vmin = vmin,vmax = vmax)
            plt.cla()

        if plt_contours:
            levels = np.arange(0.2, 3.01, 0.2)  # DON'T PLOT THE ZERO lEVEL
            CS = ax.contour(data, levels,origin='lower',linewidths=1,
                    extent=(self.x_extent[0],self.x_extent[1],
                            self.y_extent[0],self.y_extent[1]),
                    cmap = cm.bone)

            #This is not working in python 3
            # ax.clabel(CS, levels[1::2],  # label every second level
            #            inline=1,
            #            fmt='%1.1f',
            #            fontsize=10)

        # add a colorbar for the image.
        if plt_colorbar:
            cbi = plt.colorbar(im, orientation='horizontal',ax=ax,label=units)
            cbi.set_ticks(cbticks)
            cbi.ax.tick_params(axis='x')
        if title is not None:
            ax.set_title(title)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if xlabel is not None:
            ax.set_ylabel(ylabel)
        if panel_label is not None:
            ax.text(x=0.05,y=0.92,transform=ax.transAxes,s=panel_label)
        if outfile is not None:
            try:
                if outfile.endswith('png'):
                    plt.savefig(outfile, bbox_inches='tight',dpi=150)
            except:
                pass

        return fig,ax,im
        
