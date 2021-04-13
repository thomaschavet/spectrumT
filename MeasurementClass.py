from readspc import spectrum2D
import xlrd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import spectrum
from matplotlib import cm
from RandVarClass import randVar

class measurement:
    def __init__(self):#maybe put name of file in input
        sheet = xlrd.open_workbook('spectral_radiance_Wi_17G_73161PTB20.xlsx')
        sheet = sheet.sheet_by_index(0)
        L = []
        for i in range(7,sheet.nrows):
            row = [sheet.cell_value(i, 0), sheet.cell_value(i, 1), sheet.cell_value(i, 2)]
            L.append(row)
        
        frame = spectrum2D('FS_x2_500-710/FS_x2_  33-SG-700-910-Step-1-raw-Frame-01.spc')
        self.dtmeas = frame.gate_time
        self.wavelength = frame.wavelength
        len_y = len(frame.intensity)
        len_wl = len(frame.intensity[0])
        frame = spectrum2D('FS_x2_500-710_calib/gain4-SG-700-910-Step-1-raw-Frame-01.spc')
        self.dtcalib = frame.gate_time
        
        #Open file and store resuts in "intensity"
        intensity = []
        for i in range(1,10):
            if i < 10:
                num = '0' + str(i)
            else:
                num = str(i)
            frame = spectrum2D('FS_x2_500-710/FS_x2_  33-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            intensity.append(frame.intensity)
        BGintensity = []
        for i in range(1,5):
            num = str(i)
            frame = spectrum2D('FS_x2_500-710/FS_x2_bg_  33-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            BGintensity.append(frame.intensity)
        calibintensity = []
        for i in range(1,10):
            if i < 10:
                num = '0' + str(i)
            else:
                num = str(i)
            frame = spectrum2D('FS_x2_500-710_calib/gain4-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            calibintensity.append(frame.intensity)
        calibBGintensity = []
        for i in range(1,5):
            num = str(i)
            frame = spectrum2D('FS_x2_500-710_calib/gain4_bg-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            calibBGintensity.append(frame.intensity)
        '''
        pixels = []
        for j in range(9):
            pixel  = []
            for i in range(9):
                pixel.append(intensity[i][j][0].value)
            pixels.append(pixel)
        plt.plot(pixels)
        plt.show()
        '''
        #computes mean values
        self.meanU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(intensity)):
            self.meanU = self.meanU + intensity[i]
        self.meanU = self.meanU/len(intensity)
        self.meanBG = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(BGintensity)):
            self.meanBG = self.meanBG + BGintensity[i]
        self.meanBG = self.meanBG/len(BGintensity)
        self.calibmeanU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(calibintensity)):
            self.calibmeanU = self.calibmeanU + calibintensity[i]
        self.calibmeanU = self.calibmeanU/len(calibintensity)
        self.calibmeanBG = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(calibBGintensity)):
            self.calibmeanBG = self.calibmeanBG + calibBGintensity[i]
        self.calibmeanBG = self.calibmeanBG/len(calibBGintensity)
        
        #computes variances
        varU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(intensity)):
            varU = varU + (intensity[i]-self.meanU)**2
        varU = varU/len(intensity)
        varU = varU**(1/2)
        varBG = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(BGintensity)):
            varBG = varBG + (BGintensity[i]-self.meanBG)**2
        varBG = varBG/len(BGintensity)
        varBG = varBG**(1/2)
        calibvarU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(calibintensity)):
            calibvarU = calibvarU + (calibintensity[i]-self.calibmeanU)**2
        calibvarU = calibvarU/len(calibintensity)
        calibvarU = calibvarU**(1/2)
        calibvarBG = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(calibBGintensity)):
            calibvarBG = calibvarBG + (calibBGintensity[i]-self.calibmeanBG)**2
        calibvarBG = calibvarBG/len(calibBGintensity)
        calibvarBG = calibvarBG**(1/2)
        
        #Interpolate L from file for each wavelength
        self.Lcalib = []
        for wl in self.wavelength:
            Ldown = L[0][0]
            i = 0
            while Ldown < wl.value:
                 i = i + 1
                 Ldown = L[i][0]
            Lup = L[i-1][0]
            ratio = (wl.value-Ldown)/(Lup-Ldown)
            meanL = (L[i][1]*(1-ratio) + L[i-1][1]*ratio) * u.W/(u.m**3*u.sr)#u.kg/(u.m*u.s**3)
            meanL = meanL.to(u.W/(u.m**2*u.nm*u.sr)).value
            errL = L[i][2]*(1-ratio) + L[i-1][2]*ratio
            self.Lcalib.append(randVar(meanL,errL))
    
    def calibration(self,Lbound):
        Smeas = self.meanU-self.meanBG
        Scalib = self.calibmeanU-self.calibmeanBG
        
        #computes Scalib that is the signal directly measured by the plasma lamp (the lamp is at lines 530 to 540)
        meanScalib = [0 for i in range(len(Scalib[530]))]
        for i in range(530,540):
            meanScalib = meanScalib + Scalib[i]
        meanScalib = meanScalib/10
        
        #computes a random Lcalib from data : mean value and err
        randLcalib = []
        for i in range(len(self.Lcalib)):
            if Lbound == 'mean':
                randLcalib.append(self.Lcalib[i].mean())
            if Lbound == 'up':
                randLcalib.append(self.Lcalib[i].max)
            if Lbound == 'down':
                randLcalib.append(self.Lcalib[i].min)
        Lmeas = (Smeas*(1/self.dtmeas))*(meanScalib*(1/self.dtcalib))**-1 *randLcalib
        return Lmeas
    
    def spectralIntegration(self, intensity, lambda_a, lambda_b, baseline_subtraction = False):
        
        #if not isinstance(lambda_a, u.Quantity):
        #    IO.error(traceback + '"lambda_a" must be an instance of u.Quantity')
        #if not isinstance(lambda_b,u.Quantity):
        #    IO.error(traceback + '"lambda_b" must be an instance of u.Quantity')
        
        index_a = (np.abs(self.wavelength - lambda_a)).argmin()
        index_b = (np.abs(self.wavelength - lambda_b)).argmin()
        
        intensity_val = np.zeros(len(intensity))
                                      # unit = self.intensity.unit * self.wavelength.unit, 
                                      # name = 'Intensity')
        
        for i in range(len(intensity)):
            intensity_val[i] = np.trapz(x = self.wavelength[index_a:index_b+1].value, 
                                   y = intensity[i][index_a:index_b+1].value)
            
            
        intensity = intensity_val * intensity.unit * self.wavelength.unit
        dy_dp    = (170 / (922-42))  * u.mm / u.pix
        dp_dy    = 1 / dy_dp
        Dy       = 1024 * u.pix  / dp_dy #vertical dimension of probed region
        y        = np.linspace(512, -511, 1024) / 1024 * Dy
        profile = spectrum.spatialIntensity(r = y, intensity = intensity)
        
        return profile
    
    def plot2D(self, intensity):
        
        x_plt = self.wavelength

        n = intensity.shape[0]
        y_plt = np.linspace(0, n - 1, n) * u.pix
        z_plt = intensity
        
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
