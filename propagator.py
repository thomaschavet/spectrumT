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

meanP = 10000 #[pa]
errP = 0.2 #error on the pressure in %
P = randVar(meanP, errP)

lambda_a = 776.0#776.5
lambda_b = 779.0#778.0

OIlines = atomLines('O', Eu, El, meanAu, g)

#Open all image and and computes mean, variance,...
meas700_910 = measurement()

#computes mean and variance of gaussian BG
BGvect = []
for i in range(len(meas700_910.meanBG)):
    for j in range(len(meas700_910.meanBG[i])):
        BGvect.append(meas700_910.meanBG[i][j].value)
(mu, sigma) = norm.fit(BGvect)

#just to compute r and x and mean T
L700_910 = meas700_910.calibration('mean',0,0)
profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
r = profile_inverted.r
mypts = [0,2.5,5,7.5,10,15,20,25,30,40]
#mypts = [0]
x = []
for p in mypts:
    x.append(min(range(len(r)), key=lambda i: abs(r[i].value-p)))
radialT = [[] for i in range(len(x))]
meanT = []
for j in range(len(x)):
    e = emeas[x[j]]
    T = Tcalc(e.value,OIlines,P.mean())
    meanT.append(float(T))
OIlines = atomLines('O', Eu, El, Au, g)
Tmat = []
print('init done')

nTest = 1000
for i in range(nTest):
    nframe = int(random.uniform(1,11))
    ncalibframe = int(random.uniform(1,6))
    L700_910 = meas700_910.calibration('rand', nframe, ncalibframe, mu, sigma)
    profile = meas700_910.spectralIntegration(L700_910,lambda_a*u.nm,lambda_b*u.nm)
    abel_obj = abel1D.AbelInversion1D(profile)
    profile_inverted = abel_obj.invert()
    emeas = profile_inverted.intensity
    OIlines.randAu()#Take a new value for Au at each iteration
    Tvect = []
    for j in range(len(x)):
        e = emeas[x[j]]
        T = Tcalc(e.value,OIlines,P.rand())
        #T = e.value
        radialT[j].append(float(T))
        Tvect.append(float(T))
    Tmat.append(Tvect)
    if i%10==0:
        print(i)
print('computation done')

maxT = []
minT = []
for i in range(len(x)):
    maxT.append(max(radialT[i]))
    minT.append(min(radialT[i]))
#for i in range(len(Tmat)):
#    plt.fill_between(r[x],Tmat[i],meanT)#,label='O (777 nm)')
plt.fill_between(r[x],maxT,minT)
#plt.plot(r[x],radialT)
plt.xlabel('radial coordinate, mm')
plt.ylabel('Emission, W/m³')
plt.show()

plt.hist(radialT[0],bins=20,density=True)
plt.xlabel('Emission, W/m³')
plt.ylabel('Probability distribution, []')
plt.show()

plt.hist(radialT[9],bins=20,density=True)
plt.xlabel('Emission, W/m³')
plt.ylabel('Probability distribution, []')
plt.show()
