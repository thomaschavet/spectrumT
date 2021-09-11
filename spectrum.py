# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 09:16:32 2021

@author: afagn
"""

import spc
import re
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm
import IO_mod as IO

#%%=============================================================================
class spatialIntensity:
    """
    defines 1D intensity profile either from coordinate, intensity or from file
        
    required inputs
        1 arg:  
            fname -> str
        2 args: 
            wavelength  -> df.dimensionedField
            intensity   -> df.dimensionedField
        
    """
    def __init__(self, r, intensity):

        self.__declare()
        
        self.r = r
        self.intensity = intensity


    # --------------------------------------------------------------------------
    def __declare(self):
        
        self.class_name         = self.__class__.__name__

        
        self.fname = None
              
        self.r          = None    
        self.x_coord    = None
        self.y_coord    = None

        self.wavelength = None
        self.intensity  = None

        self.wavelength = None
        
        self.range      = []
        


    #---------------------------------------------------------------------------     
    def getValRadialPoint(self, R):
        
        traceback = self.class_name + '.getValRadialPoint() \n'
        
        if not isinstance(R, u.Quantity):
            IO.error(traceback + '"lambda_val" must be an instance of df.dimensionedField')
            
        intensity_point = np.interp(R.value,  
                                    self.r.value, 
                                    self.intensity.value)
        
        out = intensity_point * self.intensity.unit
        return out
            
    #---------------------------------------------------------------------------     
    def plot(self, color = 'k',
             multiplot = False, 
             marker = '.', 
             label = None):
        
        x_plt = self.r
        y_plt = self.intensity

        if not multiplot:
            plt.figure()

        plt.plot(x_plt.value, y_plt.value, color = color, marker = marker , linestyle = '--', label = label)
        plt.xlabel('Radial coordinate, '+ x_plt.unit.to_string())
        plt.ylabel('Intensity, '+ y_plt.unit.to_string())
        plt.legend()
        if not multiplot:
            plt.show() 


#%%
class spectrum1D:
    
    def __init__(self, wavelength, intensity):
        
        self.__declare()

    #--------------------------------------------------------------------------- 
    def __declare(self):   
        
        self.class_name = self.__class__.__name__

        self.wavelength = None
        self.intensity  = None
        self.gain       = None
        self.gate_time  = None
        self.accumulations = None
        
        self.y_coord    = None
        self.resolution = None
        self.range      = None 
        
    #--------------------------------------------------------------------------- 
    def baselineSubtraction(self, index_a, index_b, order = 1):

        traceback = self.class_name + '.baselineSubtraction() \n'
                
        if order == 1:
        
            n_points = 2
            
            x = np.concatenate((self.wavelength[index_a : index_a + n_points].value, 
                                self.wavelength[index_b-n_points+1 : index_b+1].value))
            y = np.concatenate((self.intensity[index_a : index_a + n_points].value, 
                                self.intensity[index_b-n_points+1 : index_b+1].value))
                
            p = np.polyfit( x, y, 1 )
            
            baseline_val = np.polyval(p, self.wavelength[index_a:index_b+1].value)
            
            baseline = baseline_val * self.intensity.unit
            
            out = self.intensity[index_a:index_b+1] - baseline
        
            
            
        else:
            IO.error(traceback + '"order" must be 1')

        return out
#%%
class spectrum2D:
    
    def __init__(self, fname):
        
        self.__declare(fname)
        self.__load()
        
        
    #--------------------------------------------------------------------------- 
    def __declare(self, fname):   
        
        self.class_name = self.__class__.__name__
        self.fname      = fname
        
        self.wavelength = None
        self.intensity  = None
        self.gain       = None
        self.gate_time  = None
        self.accumulations = None
        
        self.y_coord    = None
        self.resolution = None
        self.range      = None
        
    #--------------------------------------------------------------------------- 
    def __load(self):   
    
        fileData = spc.File(self.fname)
        
        self.__readLogContent(fileData.log_content)
        
        m = len(fileData.x)
        n = len(fileData.sub)
        
        self.wavelength         = fileData.x * u.nm
                                                     
        
        intensity_temp          = np.zeros((n, m))
                
        for i in range(n):
            intensity_temp[i, :] = fileData.sub[i].y
            
        self.intensity = intensity_temp * u.ct 
        
        self.resolution = self.intensity.shape
        self.range    = [min(self.wavelength), max(self.wavelength)]
        self.y_coord    = np.arange(0, self.resolution[0]) * u.pix
        
    #--------------------------------------------------------------------------- 
    def __readLogContent(self, log_content): 
        """
        Parameters
        ----------
        log_content : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        """
        # for i in range(len(log_content)):
        #     print(log_content[i])
      
        #-----------------------------------------------------------------------
        gain_str            = str(log_content[68])
        if 'Gain' in gain_str:
            self.gain       = int(re.findall(r'\d+', gain_str)[0])
        else: 
            print('ERROR: failed to read "Gain" from spc file log_content')
    
        #-----------------------------------------------------------------------
        repetitive_gate_str = str(log_content[74])
        if 'RepetitiveGate' in repetitive_gate_str:   
            value = float(re.findall(r'\d+', repetitive_gate_str)[1])
            self.gate_time  = (value * u.ns).to('ms')

        else:
            print('ERROR: failed to read "RepetitiveGate" from spc file log_content')
        #-----------------------------------------------------------------------
        accumu_str            = str(log_content[44])
        if 'Accumulations' in accumu_str:
            self.accumulations   = int(re.findall(r'\d+', accumu_str)[0])
        else: 
            print('ERROR: failed to read "accumulations" from spc file log_content')
    
    #--------------------------------------------------------------------------- 
    def get_lambda_index(self, lambda_val):
        
        traceback = self.class_name + '.get_lambda_index() \n'
        
        if not isinstance(lambda_val, u.Quantity):
            IO.error(traceback + '"lambda_val" must be an instance of df.dimensionedField')
            
        if lambda_val < self.wavelength[0] or lambda_val > self.wavelength[-1]:
            IO.error('lambda_val is out of spectral range')
        
        index = (np.abs(self.wavelength - lambda_val)).argmin()
        
        return index
    

    #---------------------------------------------------------------------------    
    def spectralIntegration(self, lambda_a, 
                                  lambda_b, 
                                  baseline_subtraction = False):
        
        traceback = self.class_name + '.spectralIntegration() \n'
        
        if not isinstance(lambda_a, u.Quantity):
            IO.error(traceback + '"lambda_a" must be an instance of u.Quantity')
        if not isinstance(lambda_b,u.Quantity):
            IO.error(traceback + '"lambda_b" must be an instance of u.Quantity')
            
        index_a = self.get_lambda_index(lambda_a)
        index_b = self.get_lambda_index(lambda_b)
        
        intensity_val = np.zeros(self.resolution[0])
                                      # unit = self.intensity.unit * self.wavelength.unit, 
                                      # name = 'Intensity')
                                      
        wavelength_ = self.wavelength[index_a:index_b+1]                       
        for i in range(self.resolution[0]):
            intensity_ = self.intensity[i, index_a:index_b+1]
            if baseline_subtraction:
                intensity_ = baselineSubtraction(wavelength_, 
                                                 intensity_)
            intensity_val[i] = np.trapz(x = wavelength_.value, 
                                        y = intensity_.value)
            
            
        intensity = intensity_val * self.intensity.unit * self.wavelength.unit
        
        profile = spatialIntensity(r = self.y_coord, 
                                   intensity = intensity)
            
        return profile
    
    #---------------------------------------------------------------------------      
    def plot2D(self):
        
        x_plt = self.wavelength

        n = self.intensity.shape[0]
        y_plt = np.linspace(0, n - 1, n) * u.pix
        z_plt = self.intensity
        
        xx, yy = np.meshgrid(x_plt.value, y_plt.value)
        
#        title_label = 'x = %.1f mm' % (self.x_coord)
        title_label = ''
        
        extent = [0 , 1, 
                  0 , 1]
        
        extent = [x_plt.value[0], x_plt.value[-1], 
                  y_plt.value[-1], y_plt.value[0]]


        plt.figure()
            

        
        plt.imshow(z_plt.value, 
                    cmap = cm.nipy_spectral,  
                    aspect= 'auto', 
                    extent = extent,
                    # norm=LogNorm()
                    )
        cbar = plt.colorbar()
        cbar.set_label('Intensity, ' + z_plt.unit.to_string())

        plt.xlabel('Wavelength, ' +self.wavelength.unit.to_string())
        plt.ylabel('y' +  ', ' + y_plt.unit.to_string())
        plt.title(title_label)
      
        plt.show()
            
#%%

#--------------------------------------------------------------------------- 
def baselineSubtraction(wavelength, intensity, n_points=1):
        
    

    
    x = np.concatenate((wavelength[0 : n_points].value, 
                        wavelength[-1-n_points +1 :].value))
    y = np.concatenate((intensity[0 : n_points].value, 
                        intensity[-1-n_points+1 :].value))
        
    p = np.polyfit( x, y, 1 )
    
    baseline_val = np.polyval(p, wavelength.value)
     
    baseline = baseline_val * intensity.unit
    
    out = intensity - baseline
    

    return out
