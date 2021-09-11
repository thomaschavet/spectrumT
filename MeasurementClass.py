#Open and store alle measurements. Can perform calibration of integration of one line
from readspc import spectrum2D
import xlrd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import spectrum
from matplotlib import cm
from RandVarClass import randVar
import random

class measurement:
    def __init__(self):
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
        self.intensity = []
        for i in range(1,11):
            if i < 10:
                num = '0' + str(i)
            else:
                num = str(i)
            frame = spectrum2D('FS_x2_500-710/FS_x2_  33-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            self.intensity.append(frame.intensity)
        BGintensity = []
        for i in range(1,6):
            num = str(i)
            frame = spectrum2D('FS_x2_500-710/FS_x2_bg_  33-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            BGintensity.append(frame.intensity)
        self.calibintensity = []
        for i in range(1,11):
            if i < 10:
                num = '0' + str(i)
            else:
                num = str(i)
            frame = spectrum2D('FS_x2_500-710_calib/gain4-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            self.calibintensity.append(frame.intensity)
        calibBGintensity = []
        for i in range(1,6):
            num = str(i)
            frame = spectrum2D('FS_x2_500-710_calib/gain4_bg-SG-700-910-Step-1-raw-Frame-'+num+'.spc')
            calibBGintensity.append(frame.intensity)
        
        #computes mean values
        self.meanU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(self.intensity)):
            self.meanU = self.meanU + self.intensity[i]
        self.meanU = self.meanU/len(self.intensity)
        self.meanBG = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(BGintensity)):
            self.meanBG = self.meanBG + BGintensity[i]
        self.meanBG = self.meanBG/len(BGintensity)
        self.calibmeanU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(self.calibintensity)):
            self.calibmeanU = self.calibmeanU + self.calibintensity[i]
        self.calibmeanU = self.calibmeanU/len(self.calibintensity)
        self.calibmeanBG = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(calibBGintensity)):
            self.calibmeanBG = self.calibmeanBG + calibBGintensity[i]
        self.calibmeanBG = self.calibmeanBG/len(calibBGintensity)
        
        #computes standard deviation
        self.stdU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(self.intensity)):
            self.stdU = self.stdU + (self.intensity[i]-self.meanU)**2
        self.stdU = self.stdU/len(self.intensity)
        self.stdU = self.stdU**(1/2)
        self.calibstdU = [[0 for col in range(len_wl)] for row in range(len_y)]
        for i in range(len(self.calibintensity)):
            self.calibstdU = self.calibstdU + (self.calibintensity[i]-self.calibmeanU)**2
        self.calibstdU = self.calibstdU/len(self.calibintensity)
        self.calibstdU = self.calibstdU**(1/2)
        
        #Interpolate L from file for each wavelength
        self.Lcalib = []
        self.vectmeanL = []
        self.stdDev = []
        for j in range(len(self.wavelength)):
            Ldown = L[0][0]
            i = 0
            while Ldown < self.wavelength[j].value:
                 i = i + 1
                 Ldown = L[i][0]
            Lup = L[i-1][0]
            ratio = (self.wavelength[j].value-Ldown)/(Lup-Ldown)
            meanL = (L[i][1]*(1-ratio) + L[i-1][1]*ratio) * u.W/(u.m**3*u.sr)
            meanL = meanL.to(u.W/(u.m**2*u.nm*u.sr)).value
            errL = L[i][2]*(1-ratio) + L[i-1][2]*ratio
            self.vectmeanL.append(meanL)
            self.stdDev.append(errL/2) #err is 2*standard deviation
            self.Lcalib.append(randVar(meanL,errL))
        stdDev = np.array(self.stdDev)
        self.covL = np.outer(stdDev,stdDev.T)
    
    def calibration(self,Lbound, nframe, ncalibframe, mu='none', sigma='none'):
        #nframe is the num of the frame that should be used for the computation
        #if nframe == 0 we take the mean frame
        #if nframe < 0 we generate a frame from mean and standard deviation
        if nframe == 0:
            image = self.meanU
        elif nframe < 0:
            corr = 'correlated'
            if corr == 'uncorrelated':
                stddev = np.random.normal(0, 1, [len(self.meanU),len(self.meanU[0])])
                image = self.meanU + np.multiply(stddev,self.stdU)
            elif corr == 'correlated':
                stddev = np.random.normal(0, 1)
                image = self.meanU + self.stdU*stddev
        else:
            image = self.intensity[nframe-1]
        if ncalibframe == 0:
            calibimage = self.calibmeanU
        elif ncalibframe < 0:
            corr = 'correlated'
            if corr == 'uncorrelated':
                stddev = np.random.normal(0, 1, [len(self.calibmeanU),len(self.calibmeanU[0])])
                calibimage = self.calibmeanU + np.multiply(stddev,self.calibstdU)
            elif corr == 'correlated':
                stddev = np.random.normal(0, 1)
                calibimage = self.calibmeanU + self.calibstdU*stddev
        else:
            calibimage = self.calibintensity[nframe-1]
        
        #subtraction of the BG image
        if mu == 'none' and sigma == 'none':
            Smeas = image-self.meanBG
            Scalib = calibimage-self.calibmeanBG
        else:
            genBG = np.random.normal(mu, sigma, [len(self.meanBG),len(self.meanBG[0])]) *u.ct
            Smeas = image-genBG
            genBG = np.random.normal(mu, sigma, [len(self.calibmeanBG),len(self.calibmeanBG[0])]) *u.ct
            Scalib = calibimage-genBG
        
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
            if Lbound == 'rand':
                randLcalib.append(self.Lcalib[i].rand())
        if Lbound == 'randcorr':#random correleted
            var = np.random.normal()
            randLcalib = self.vectmeanL + np.array(self.stdDev)*var
        self.Lmeas = (Smeas*(1/self.dtmeas))*(meanScalib*(1/self.dtcalib))**-1 *randLcalib
        return self.Lmeas
    
    def newlambda(self, lambda_a, lambda_b,scale):
        #computes the two integration bound that correspond to lambda at maximum of the peak + or - scale*with_at_half_max
        index_a = (np.abs(self.wavelength.value - lambda_a)).argmin()
        index_b = (np.abs(self.wavelength.value - lambda_b)).argmin()
        centerline = int(len(self.Lmeas)/2)
        highpeak = max(self.Lmeas[centerline][index_a:index_b])
        indexpeak = (np.abs(self.Lmeas[centerline] - highpeak)).argmin()
        centerpeak = self.wavelength.value[indexpeak]
        halfmax = highpeak/2
        
        i = indexpeak
        while halfmax < self.Lmeas[centerline][i]:
            ratio = (halfmax - self.Lmeas[centerline][i-1])/(self.Lmeas[centerline][i]-self.Lmeas[centerline][i-1])
            i = i - 1
        wl1 = self.wavelength.value[i] + ratio * (self.wavelength.value[i+1]-self.wavelength.value[i])
        i = indexpeak
        while halfmax < self.Lmeas[centerline][i]:
            ratio = (self.Lmeas[centerline][i]-halfmax)/(self.Lmeas[centerline][i]-self.Lmeas[centerline][i+1])
            i = i + 1
        wl2 = self.wavelength.value[i-1] + ratio * (self.wavelength.value[i]-self.wavelength.value[i-1])
        HMW = wl2 - wl1
        lambda_a = centerpeak - scale*HMW
        lambda_b = centerpeak + scale*HMW
        return lambda_a,lambda_b
    
    
    def spectralIntegration(self, intensity, lambda_a, lambda_b, errSpatialCalib = False):
        #integration of one emission line and spatial calibration
        index_a = (np.abs(self.wavelength - lambda_a)).argmin()
        index_b = (np.abs(self.wavelength - lambda_b)).argmin()
        
        intensity_val = np.zeros(len(intensity))
        
        wavelength_ = self.wavelength[index_a:index_b+1]
        for i in range(len(intensity)):
            intensity_ = intensity[i][index_a:index_b+1]
            intensity_ = self.baselineSubtraction(wavelength_, intensity_)#subtract baseline
            intensity_val[i] = np.trapz(x = wavelength_.value, 
                                   y = intensity_.value)
        
        intensity = intensity_val * intensity.unit * self.wavelength.unit
        res=1#pixel
        if errSpatialCalib == True:
            dy_dp    = (170 / (922+random.uniform(-res/2, res/2)-42+random.uniform(-res/2,res/2)))  * u.mm / u.pix
        else:
             dy_dp    = (170 / (922-42))  * u.mm / u.pix
        dp_dy    = 1 / dy_dp
        Dy       = 1024 * u.pix  / dp_dy #vertical dimension of probed region
        leny = len(intensity)
        y        = np.linspace(int(leny/2), -int(leny/2), leny) / leny * Dy
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
        cbar.set_label('Number of photons, []')

        plt.xlabel('Wavelength, [' +self.wavelength.unit.to_string() + ']')
        plt.ylabel('y, [' + y_plt.unit.to_string() + ']')
        plt.title(title_label)
        plt.savefig('Terr.pdf', format='pdf')
        plt.show()
        
    def baselineSubtraction(self, wavelength, intensity, n_points=5):
        
        x = np.concatenate((wavelength[0 : n_points].value, 
                            wavelength[-1-n_points +1 :].value))
        y = np.concatenate((intensity[0 : n_points].value, 
                            intensity[-1-n_points+1 :].value))
            
        p = np.polyfit( x, y, 1 )
        
        baseline_val = np.polyval(p, wavelength.value)
         
        baseline = baseline_val * intensity.unit
        
        out = intensity - baseline
        
    
        return out