from RandVarClass import randVar
from AtomicLinesClass import atomLines
from function import efromSpectr
from function import Tcalc
import matplotlib.pyplot as plt
from MeasurementClass import measurement
import astropy.units as u
#import spectrum
import abel1D

#Definition of the atomic lines (Here OI between 776 and 778 nm)
mtoJ = 1.986445857E-25
Eu = [8663145.4*mtoJ,8662777.8*mtoJ,8662575.7*mtoJ] #from [m^-1] to [J]
El = 7376820*mtoJ
g = [7,5,3]
meanAu = 3.69E7
Au = randVar(meanAu,0)#3)
OIlines = atomLines('O', Eu, El, Au, g)

meanP = 10000 #[pa]
errP = 0#.2 #error on the pressure in %
P = randVar(meanP, errP)

lambda_a = 776.5
lambda_b = 778.0

#emeas = efromSpectr(lambda_a, lambda_b)
#print(emeas)

meas700_910 = measurement()
L700_910 = meas700_910.calibration()
#meas700_910.plot2D(L700_910)
#plt.plot(L700_910[512])
#plt.show()
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
#profile.plot()
abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
#print(emeas[0].value)

#profile_inverted.plot()
#x = []
#for i in range(len(L700_910)):
#    x.append(L700_910[i][0])
#plt.plot(x)
#plt.show()
#plt.imshow(L700_910)
#plt.show()
#profile = spectrum.spatialIntensity(r = y_coord, intensity = intensity)
#abel_obj = abel1D.AbelInversion1D(profile)
#profile_inverted = abel_obj.invert()
'''
BGvect = []
for i in range(len(meas700_910.meanBG)):
    for j in range(len(meas700_910.meanBG[i])):
        if 400 < meas700_910.meanBG[i][j].value < 700:
            BGvect.append(meas700_910.meanBG[i][j].value)
plt.hist(BGvect)
plt.show()
'''
e = randVar(emeas[0].value, 0)

results = []
nTest = 1
for i in range(nTest):
    OIlines.randAu()#Take a new value for Au at each iteration
    T = Tcalc(e.rand(),OIlines,P.rand())
    results.append(float(T))
    if i%100==0:
        print(i)
print(results)
#plt.hist(results,density=True)
#plt.xlabel('Temperature')
#plt.ylabel('probability distribution')
#plt.show()
