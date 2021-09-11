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

#%%
class spectrum2D:
    
    def __init__(self, fname):
        
        self.wavelength = None
        self.intensity  = None
        self.gain       = None
        self.gate_time  = None
        
        
        fileData = spc.File(fname)
        
        self.__readLogContent(fileData.log_content)
        
        m = len(fileData.x)
        n = len(fileData.sub)
        
        self.wavelength         = fileData.x * u.nm
                                                     
        
        intensity_temp          = np.zeros((n, m))
                
        for i in range(n):
            intensity_temp[i, :] = fileData.sub[i].y
            
        self.intensity = intensity_temp * u.ct 
        
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
#fname = 'D:/VKI/PhD/jobs/SPECTRO/mirrors/2021_01_28/HF2000/spc/FS_x2_  33-SG-700-910-Step-1-raw-Frame-01.spc'
#fname = '/home/thomas/Desktop/FS_x2_  33-SG-700-910-Step-1-raw-Frame-01.spc'

#test = spectrum2D(fname)

#test.plot2D()
#print(test.gain)
#print(test.gate_time)
