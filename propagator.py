from RandVarClass import randVar
from AtomicLinesClass import atomLines
from function import Tcalc
import matplotlib.pyplot as plt
from MeasurementClass import measurement
import astropy.units as u
import abel1D
from scipy.stats import norm
import random

#Definition of the atomic lines (Here OI between 776 and 778 nm) Data from NIST
mtoJ = 1.986445857E-25
Eu = [8663145.4*mtoJ,8662777.8*mtoJ,8662575.7*mtoJ] #from [m^-1] to [J]
El = 7376820*mtoJ
g = [7,5,3]
meanAu = 3.69E7
Au = randVar(meanAu,3)
meanAu = randVar(Au.mean(),0)

meanP = 10000 #[pa]
errP = 2 #error on the pressure in %
P = randVar(meanP, errP)

lambda_a = 774#nm
lambda_b = 782#nm

OIlines = atomLines('O', Eu, El, meanAu, g)

#Open all image and computes mean, variance,...
meas700_910 = measurement()

#computes mean and variance of gaussian BG
BGvect = []
for i in range(len(meas700_910.meanBG)):
    for j in range(len(meas700_910.meanBG[i])):
        BGvect.append(meas700_910.meanBG[i][j].value)
(mu, sigma) = norm.fit(BGvect)

#computes the mean temperature
L700_910 = meas700_910.calibration('mean',0,0)#calibration
[lambda_a,lambda_b] = meas700_910.newlambda(lambda_a,lambda_b,3)#computes bound of integration
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)#integrate emission line
abel_obj = abel1D.AbelInversion1D(profile)#initiate Abel inversion
profile_inverted = abel_obj.invert()#Abel inversion
emeas = profile_inverted.intensity
r = profile_inverted.r.value

#define at which radial coordinate do we compute the temperature
mypts = [0,2.5,5,7.5,10,15,20,25,30,40,50]#mm
x = []
for p in mypts:
    x.append(min(range(len(r)), key=lambda i: abs(r[i]-p)))

#computes T at different radial coordinates
radialTO = [[] for i in range(len(x))]
meanTO = []
for j in range(len(x)):
    e = emeas[x[j]]
    T = Tcalc(e.value,OIlines,P.mean())
    meanTO.append(float(T))


#what parameters to take into account
#put True to consider it and Fasle to ignore it
errAu = True
errP = True
errImage = True
errCalibImage = True
errRadCalib = True
errBG = True
errBL = True
errSpatialCalib = True

print('init done')
#here we recompute the temperature ntest number of time by taking uncertainties into account

if errAu == True:
    OIlines = atomLines('O', Eu, El, Au, g)
else:
    OIlines = atomLines('O', Eu, El, meanAu, g)
nTest = 1000
for i in range(nTest):
    if errBL == True:
        newlambda_a = lambda_a+random.uniform(-0.3, 0.3)
        newlambda_b = lambda_b+random.uniform(-0.3, 0.3)
    else:
        newlambda_a = lambda_a
        newlambda_b = lambda_b
    
    if errImage == True:
        nframe = -1# it means we generate a random image from mean and std deviation
    else:
        nframe = 0#it means we take the mean image
    if errCalibImage == True:
        ncalibframe = -1# it means we generate a random image from mean and std deviation
    else:
        ncalibframe = 0#it means we take the mean image
    
    if errRadCalib == True:
        radcalib = 'randcorr'#put 'rand' for uncorrelated and 'randcorr' for fully correlated
    else:
        radcalib = 'mean'
    
    if errP == True:
        myP = P.rand()
    else:
        myP = P.mean()
    
    if errBG == True:
        L700_910 = meas700_910.calibration(radcalib, nframe, ncalibframe, mu, sigma)
    else:
        L700_910 = meas700_910.calibration(radcalib, nframe, ncalibframe)
    
    profile = meas700_910.spectralIntegration(L700_910,newlambda_a*u.nm,newlambda_b*u.nm,errSpatialCalib)
    abel_obj = abel1D.AbelInversion1D(profile)
    profile_inverted = abel_obj.invert()
    emeas = profile_inverted.intensity
    
    OIlines.randAu()#Take a new value for Au at each iteration
    Tvect = []
    for j in range(len(x)):
        e = emeas[x[j]].value
        if e<0:
            e=0
        T = Tcalc(e,OIlines,myP)
        radialTO[j].append(float(T))
        Tvect.append(float(T))
    
    if i%10==0:
        print(i)
print('computation done')
maxTO = []
minTO = []
for i in range(len(x)):
    maxTO.append(max(radialTO[i]))
    minTO.append(min(radialTO[i]))


#now the exact same code is used but considering the nitrogen line at 747 nm


#Definition of the atomic lines (Here NI between 746 and 748 nm)
mtoJ = 1.986445857E-25
Eu = [9675084 * mtoJ] #from [m^-1] to [J]
El = [8336462 * mtoJ]
g = [4]
meanAu = 1.96E7
Au = randVar(meanAu,3)
meanAu = randVar(Au.mean(),0)

meanP = 10000 #[pa]
errP = 2 #error on the pressure in %
P = randVar(meanP, errP)

lambda_a = 746#746.2
lambda_b = 748.5#747.4

NIlines = atomLines('N', Eu, El, meanAu, g)

#Open all image and and computes mean, variance,...
meas700_910 = measurement()

#computes mean and variance of gaussian BG
BGvect = []
for i in range(len(meas700_910.meanBG)):
    for j in range(len(meas700_910.meanBG[i])):
        BGvect.append(meas700_910.meanBG[i][j].value)
(mu, sigma) = norm.fit(BGvect)

L700_910 = meas700_910.calibration('mean',0,0)
[lambda_a,lambda_b] = meas700_910.newlambda(lambda_a,lambda_b,3)
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
r = profile_inverted.r.value

mypts = [0,2.5,5,7.5,10,15,20,25,30,40,50]
x = []
for p in mypts:
    x.append(min(range(len(r)), key=lambda i: abs(r[i]-p)))

radialTN = [[] for i in range(len(x))]
meanTN = []
for j in range(len(x)):
    e = emeas[x[j]]
    T = Tcalc(e.value,NIlines,P.mean())
    meanTN.append(float(T))

#what parameters to take into account
errAu = True
errP = True
errImage = True
errCalibImage = True
errRadCalib = True
errBG = True
errBL = True
errSpatialCalib = True

print('init done')

if errAu == True:
    NIlines = atomLines('N', Eu, El, Au, g)
else:
    NIlines = atomLines('N', Eu, El, meanAu, g)
nTest = 1000
for i in range(nTest):
    if errBL == True:
        newlambda_a = lambda_a+random.uniform(-0.3, 0.3)
        newlambda_b = lambda_b+random.uniform(-0.3, 0.3)
    else:
        newlambda_a = lambda_a
        newlambda_b = lambda_b
    
    if errImage == True:
        nframe = -1# it means we generate a random image from mean and std deviation
    else:
        nframe = 0#it means we take the mean image
    if errCalibImage == True:
        ncalibframe = -1# it means we generate a random image from mean and std deviation
    else:
        ncalibframe = 0#it means we take the mean image
    
    if errRadCalib == True:
        radcalib = 'randcorr'#put 'rand' for uncorrelated and 'randcorr' for fully correlated
    else:
        radcalib = 'mean'
    
    if errP == True:
        myP = P.rand()
    else:
        myP = P.mean()
    
    if errBG == True:
        L700_910 = meas700_910.calibration(radcalib, nframe, ncalibframe, mu, sigma)
    else:
        L700_910 = meas700_910.calibration(radcalib, nframe, ncalibframe)
    
    profile = meas700_910.spectralIntegration(L700_910,newlambda_a*u.nm,newlambda_b*u.nm,errSpatialCalib)
    abel_obj = abel1D.AbelInversion1D(profile)
    profile_inverted = abel_obj.invert()
    emeas = profile_inverted.intensity
    
    NIlines.randAu()#Take a new value for Au at each iteration
    Tvect = []
    for j in range(len(x)):
        e = emeas[x[j]].value
        if e<0:
            e=0
        T = Tcalc(e,NIlines,myP)
        #T = e.value
        radialTN[j].append(float(T))
        Tvect.append(float(T))
    
    if i%10==0:
        print(i)
print('computation done')

maxTN = []
minTN = []
for i in range(len(x)):
    maxTN.append(max(radialTN[i]))
    minTN.append(min(radialTN[i]))

#plot the temperature and uncertainties
plt.fill_between(r[x],maxTO,minTO,color=(0.5,0.8,1),alpha=1)
plt.plot(r[x],meanTO,color=(0.2,0.5,0.7))
plt.fill_between(r[x],maxTN,minTN,color=(1,0.8,0.5),alpha=0.6)
plt.plot(r[x],meanTN,color=(0.7,0.5,0.2))
plt.xlabel('Radial coordinate, [mm]')
plt.ylabel('Temperature, [K]')
plt.legend(['Oxygen triplet 777 nm','Nitrogen line 747 nm'])
#plt.savefig('Terr.pdf', format='pdf')
plt.show()

#plot histogram of T at the center of the jet considering O line 777 nm
n, bins, patches = plt.hist(radialTO[0],bins=20,density=True)
(mu, sigma) = norm.fit(radialTO[0])
y = norm.pdf(bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)
plt.xlabel('Temperature, [K]')
plt.ylabel('Probability distribution, [K⁻¹]')
#plt.savefig('hist0O.pdf', format='pdf')
plt.show()
print('mean T: ',mu)
print('error: +-',sigma*2)