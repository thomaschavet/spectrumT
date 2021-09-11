# -*- coding: utf-8 -*-



#%%
#import spectrum2D
import abel
import astropy.units as u
import IO_mod as IO
import utilities as ut
import numpy as np
from scipy.interpolate import interp1d
import spectrum
#%% ======================================================================================================================

class butterParameters:
    def __init__(self):
        self.fc = 10
        self.n  = 6
            
class convParameters:
    def __init__(self):
        self.box_pnts = 50
         
class filterParameters:
    def __init__(self):
        self.method     = 'butter'
        self.butter     = butterParameters()
        self.conv       = convParameters()
        
class parameters:
    method              = 'hansenlaw' # 'basex'/ 'hansenlaw'
    R_max               = 90 * u.mm # max distance from determined centerline
    n_points            = 501
    n_pnts_disc          = 1
    filtering           = filterParameters()
#%% ======================================================================================================================
class AbelInversion1D:
    """
    Computes Abel inversion of the double-sided (top and bottom) intensity profile I = I(y)
    If not defined, it computes the profile center and averages the top and bottom 
    sides before inversion.
    The profile is smoothed with a butterworth filter before inversion.
    """
    parameters = parameters()
    
    # method              = 'hansenlaw' # 'basex'/ 'hansenlaw'
    # R_max               = df.dimensionedField(value = 80, unit = 'mm') # max distance from determined centerline
    # n_points            = 511
    # n_pnts_disc         = 0
    # filtering           = filterParameters()
    
    def __init__(self, spatial_intensity, **kwargs):

        self.__declare()
        
        self.y          = spatial_intensity.r
        self.intensity  = spatial_intensity.intensity
        
        self.__check()
        
        if 'y0_index' in kwargs:
            self.y0_index   = int(kwargs.get('y0_index', None))
            


    #--------------------------------------------------------------------------
    def __declare(self):

        self.debug_plot     = False

        self.y              = None # vertical coordinate
        self.intensity      = None # intensity

        self.y0_index       = None # index of the profile center

        self.r              = np.linspace(0, parameters.R_max.value, parameters.n_points) * parameters.R_max.unit #  radial coordinate
        
        self.dr             = (self.r[1] - self.r[0]) # for equally spaced

    #---------------------------------------------------------------------------      
    def __check(self):
        
        if self.y.shape != self.intensity.shape:
            IO.error('shapes of attributes "y" and "intensity" are not compatible')
            
        if not isinstance(self.y, u.Quantity):
            IO.error('"y" must be an instance of astropy.units.Quantity')
            
        if not isinstance(self.intensity, u.Quantity):
            IO.error('"intensity" must be an instance of astropy.units.Quantity')
        
    #---------------------------------------------------------------------------      
    def smooth(self, method = None):
        
        if method is None:
            method = self.parameters.filtering.method
            

        if method == 'conv':
            box_pnts = self.parameters.filtering.conv.box_pnts
            self.intensity_smooth = ut.filtering_1D(self.y.value, self.intensity.value).smooth_conv(box_pnts)
        elif method == 'butter':
            fs  = self.y.size  # sampling frequency
            fc  = self.parameters.filtering.butter.fc # cut-off frequency
            n   = self.parameters.filtering.butter.n  # filter order
            
            self.intensity_smooth = ut.filtering_1D(self.y, 
                                                    self.intensity).butter( fs, fc, n ) 
            # print(self.intensity_smooth)

        if self.debug_plot:
           plt.figure()
           plt.plot(self.y, self.intensity, label = 'original')
           plt.plot(self.y, self.intensity_smooth, label = 'filtered')
           plt.legend()
           plt.show()
        
        return self.intensity_smooth
    
    #--------------------------------------------------------------------------- 
    def find_profile_center(self):
        
        
        y_fit       = ut.fitting_1D(self.y.value, self.intensity_smooth.value).gaussian()
        y0_index    = y_fit.argmax()

        if self.debug_plot:
           plt.figure()
           plt.plot(self.y, self.intensity_smooth,'b+:',label='data')
           plt.plot(self.y, y_fit)
           plt.legend()
           plt.show()
        
        self.y0_index = y0_index
        
        return y0_index

    # ---------------------------------------------------------------------------
    def average_top_bottom(self):
        # debug_plot = True
        # print(self.n_pnts_disc)
        # print(self.intensity_smooth.shape)
        # print(self.y0_index)
        # get top and bottom half-profiles
        intensity_top = np.flipud(self.intensity_smooth[0 : self.y0_index+1])
        intensity_bottom = self.intensity_smooth[self.y0_index :]

        y_top = np.flipud(self.y[0 : self.y0_index+1])
        y_bottom = self.y[self.y0_index :]

        # eliminate n_pnts_discard from edges
        
        intensity_top_disc = intensity_top[0:- self.parameters.n_pnts_disc]
        intensity_bottom_disc = intensity_bottom[0:- self.parameters.n_pnts_disc]

        y_top_disc = y_top[0:- self.parameters.n_pnts_disc]
        y_bottom_disc = y_bottom[0:- self.parameters.n_pnts_disc]
        
        # print(intensity_top_disc.shape)
        # print(intensity_bottom_disc.shape)
        # print(y_top_disc.shape)
        # print(y_bottom_disc.shape)
        
        # print(intensity_bottom_disc[0] - intensity_top_disc[0])

        # print(y_top_disc[0])
        # print(y_bottom_disc[0])
        
        # pad shortest half-profile to compute average
        if intensity_bottom_disc.size > intensity_top_disc.size:
            n_pnts = intensity_top_disc.size
            intensity_avg = np.add(intensity_bottom_disc[0:n_pnts], intensity_top_disc) / 2.
            
            # print(intensity_bottom_disc[0:n_pnts] - intensity_top_disc)
            
            y_avg = - (y_bottom_disc[0:n_pnts] - self.y[self.y0_index])
        elif intensity_bottom_disc.size < intensity_top_disc.size:
            n_pnts = intensity_bottom_disc.size
            # print(intensity_top_disc[0:n_pnts] - intensity_bottom_disc)
            intensity_avg = np.add(intensity_top_disc[0:n_pnts], intensity_bottom_disc) / 2.
            y_avg = y_top_disc[0:n_pnts] - self.y[self.y0_index]
            # print(y_top_disc[0])
            # print(self.y[self.y0_index])
        else:
            intensity_avg = np.add(intensity_top_disc, intensity_bottom_disc) / 2.
            y_avg = y_top_disc - self.y[self.y0_index]
            
        # print(intensity_avg.shape)
        # print(y_avg)
        

        f = interp1d(y_avg, intensity_avg)
        #        intensity_avg_interp = np.interp(y_lin, y_avg, intensity_avg )
        # print(self.r)

        if self.debug_plot:
            plt.figure()
            plt.plot(y_top, intensity_top, label = 'top')
            plt.plot(-y_bottom, intensity_bottom, label = 'bottom')
            plt.plot(y_avg, intensity_avg, label = 'average')
            plt.legend()
            plt.show()

        #        self.intensity_avg  = intensity_avg_interp
        # print(y_avg)
        # print(self.r.shape)
        self.intensity_avg = f(self.r) * self.intensity.unit
        # self.intensity_avg = self.intensity_avg - self.intensity_avg[-1]
        
        return intensity_avg

    
    #---------------------------------------------------------------------------
    def impose_zero_derivative(self):     
        
        intensity_to_invert     = self.intensity_avg
        intensity_to_invert[0]  = intensity_to_invert[1]
        
        self.intensity_to_invert = intensity_to_invert

    #---------------------------------------------------------------------------
    def preprocess(self):
        
        self.smooth()
        
        if self.y0_index is None:
            self.find_profile_center()
            # print(self.y0_index)
            
        self.average_top_bottom()

#        self.impose_zero_derivative()
        
        return self.r, self.intensity_avg
        
    #---------------------------------------------------------------------------
    def invert(self, method = None):

        if method is None:
            method = self.parameters.method
             
        self.preprocess()

        dr_ = self.dr.to('m')
        
        if method == 'hansenlaw':
            intensity_inverted = abel.hansenlaw.hansenlaw_transform(self.intensity_avg.value, 
                                                                    direction='inverse', 
                                                                    dr = dr_.value)
        elif method == 'basex':
            intensity_inverted = abel.basex.basex_transform(self.intensity_avg.value, 
                                                            direction='inverse',
                                                            dr = dr_.value) 

        intensity_df = intensity_inverted * self.intensity_avg.unit / dr_.unit
        
#        out = out.convert_to('W/m3/sr/nm')
        
        profile = spectrum.spatialIntensity(r            = self.r, 
                                   intensity    = intensity_df)
        
        return profile

