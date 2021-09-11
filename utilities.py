# -*- coding: utf-8 -*-

from scipy import signal
from scipy.optimize import curve_fit
import numpy as np
import astropy.units as u
# %% --------------------------------------------------------------------------


class fitting_1D:
    
    def __init__(self, x, y):
        
        self.x = x
        self.y = y
        
        
    def gaussian(self):
        
        method = 'curve_fit'
        
        if method == 'old':
        
            y_max = self.y.max()    
            y_norm = self.y / y_max
            
            X = np.arange(self.y.size)
            
            x = np.sum(X*self.y)/np.sum(self.y)
            width = np.sqrt(np.abs(np.sum((X-x)**2*self.y)/np.sum(self.y)))
            
            max = self.y.max()
            
            fit = lambda t : max*np.exp(-(t-x)**2/(2*width**2))
            
            y_fit = fit(X)
            
            
        elif method == 'curve_fit': 
            
            if isinstance(self.x, u.Quantity):
                x_ = self.x.value
            else:
                x_ = self.x
                
            if isinstance(self.y, u.Quantity):
                y_ = self.y.value
            else:
                y_ = self.y
        
            def Gauss(x, a, x0, sigma):
                return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
            
            mean = np.sum(x_ * y_) / np.sum(y_)
            sigma = np.sqrt(abs(np.sum(y_ * (x_- mean)**2) / np.sum(y_)))

            popt, pcov = curve_fit(Gauss, x_, y_, p0=[np.max(y_), mean, sigma])

            y_fit = Gauss(x_, *popt)

        
#        fitter          = modeling.fitting.LevMarLSQFitter()
#        model           = modeling.models.Gaussian1D()   # depending on the data you need to give some initial values
#        fitted_model    = fitter(model, self.x, y_norm)
        
#        y_fit   = fitted_model(self.x) * y_max
        
        self.y_fit = y_fit
        
        return y_fit
    
# %% -------------------------------------------------------------------------


class filtering_1D:
    
    def __init__(self, x, y):
        
        self.x = x
        self.y = y
        
        if np.isnan(np.sum(self.y)):
            IO.error('y has nan values')
        
    def smooth_conv(self, box_pnts):
        """
        smooth signal by convolving with a box function (same as moving average)
        """
        
        box = np.ones(box_pnts)/box_pnts
        
        y_filt_val = np.convolve(self.y, box, mode='same')
        
        if isinstance(self.y, df.dimensionedField):
            y_filt = df.dimensionedField(value = y_filt_val, 
                                         unit = self.y.unit, 
                                         name = self.y.getName())
        else:
            y_filt = y_filt_val
        
        return y_filt
    
    def butter(self, fs, fc, n):
        """
        smooth signal with butterworth filter
        
        fs  = sampling frequency
        fc  = cut-off frequency
        n   = filter order
        """
        
        w = fc / (fs / 2) # Normalize the frequency
        
        b, a = signal.butter(n, w, 'low')
        
        y_filt_val = signal.filtfilt(b, a, self.y)
        
        if isinstance(self.y, u.Quantity):
            y_filt =  y_filt_val *self.y.unit
        else:
            y_filt = y_filt_val
            
        
        return y_filt