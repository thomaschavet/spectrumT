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
Au = randVar(meanAu,3)
minAu = randVar(Au.min,0)
maxAu = randVar(Au.max,0)
meanAu = randVar(Au.mean(),0)
OIlines = atomLines('O', Eu, El, meanAu, g)

meanP = 10000 #[pa]
errP = 0.2 #error on the pressure in %
P = randVar(meanP, errP)

lambda_a = 776.0#776.5
lambda_b = 779.0#778.0

meas700_910 = measurement()
L700_910 = meas700_910.calibration('mean')
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
r = profile_inverted.r
'''
#Plot histogram of pixel (for background measurement)
BGvect = []
for i in range(len(meas700_910.meanBG)):
    for j in range(len(meas700_910.meanBG[i])):
        if 400 < meas700_910.meanBG[i][j].value < 700:
            BGvect.append(meas700_910.meanBG[i][j].value)
n, bins, patches = plt.hist(BGvect,bins=1060,density=True)
(mu, sigma) = norm.fit(BGvect)
y = norm.pdf( bins, mu, sigma)#*len(BGvect)
l = plt.plot(bins, y, 'r--', linewidth=2)
plt.xlabel('signal measured by pixels')
plt.ylabel('density')
plt.show()
print(mu)
print(sigma)
'''
mypts = [0,2.5,5,7.5,10,15,20,25,30,35,40]
x = []
for p in mypts:
    x.append(min(range(len(r)), key=lambda i: abs(r[i].value-p)))

#"True" T
radialT = []
for j in x:
    e = emeas[j]
    results = []
    nTest = 1
    for i in range(nTest):
        OIlines.randAu()#Take a new value for Au at each iteration
        T = Tcalc(e.value,OIlines,P.mean())
        results.append(float(T))
        #if i%100==0:
            #print(i)
    radialT.append(results)
plt.plot(r[x],radialT,label='O (777 nm)')
centerT = radialT[0][0]

'''
#plot temperature computed from nitrogen line
Eu = 9675084 * mtoJ #from [m^-1] to [J]
El = 8336462 * mtoJ
g = 4
Au = randVar(1.96E7,0)
NIlines = atomLines('N', Eu, El, Au, g)
lambda_a = 746.2
lambda_b = 747.4
meas700_910 = measurement()
L700_910 = meas700_910.calibration('mean')
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
r = profile_inverted.r
radialT = []
for j in x:
    e = emeas[j]
    T = Tcalc(e.value,NIlines,P.mean())
    radialT.append(float(T))
plt.plot(r[x],radialT,label='N (747 nm)')
'''

#max T
L700_910 = meas700_910.calibration('up')
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
r = profile_inverted.r
OIlines = atomLines('O', Eu, El, minAu, g)
radialT = []
for j in x:
    e = emeas[j]
    T = Tcalc(e.value,OIlines,P.min)
    radialT.append(float(T))
plt.plot(r[x],radialT)
centerTmax = radialT[0]

#min T
L700_910 = meas700_910.calibration('down')
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
r = profile_inverted.r
OIlines = atomLines('O', Eu, El, maxAu, g)
radialT = []
for j in x:
    e = emeas[j]
    T = Tcalc(e.value,OIlines,P.max)
    radialT.append(float(T))
plt.plot(r[x],radialT)
centerTmin = radialT[0]

plt.xlabel('radial coordinate, mm')
plt.ylabel('Temperature, K')
plt.legend()
#plt.savefig('FinalTbound.eps', format='eps')
plt.show()

print(centerT)
print(centerTmax)
print(centerTmin)
print((centerTmax-centerT)/centerT*100)
print((centerTmin-centerT)/centerT*100)
