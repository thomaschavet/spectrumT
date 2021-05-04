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
plt.plot(profile.r.to(u.m),profile.intensity)
start = 100
end = 550
def poly(x,a0,a2,a3,a4,a5,a6,a7,a8):#we want a polynome with zero slope at 0 so a1=0
    return a0+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8
#a = np.polyfit(profile.r[start:end].to(u.m).value,profile.intensity[start:end].value,8)
a = curve_fit(poly, profile.r[start:end].to(u.m).value, profile.intensity[start:end].value)[0]
a = np.insert(a,1,0)#a1=0

#plt.plot(profile.r[0:end],[a[9]+a[8]*x+a[7]*x**2+a[6]*x**3+a[5]*x**4+a[4]*x**5+a[3]*x**6+a[2]*x**7+a[1]*x**8+a[0]*x**9 for x in profile.r[0:end].value])
#plt.show()
#a0 = a[-1]
#print(a[8])
r = np.flip(profile.r[200:512].to(u.m).value)
AbelMatrix = np.array(AbelInvMatrix(r,8))
plt.plot([x for x in np.arange(0,0.07,0.001)],[a[0]+a[1]*x+a[2]*x**2+a[3]*x**3+a[4]*x**4+a[5]*x**5+a[6]*x**6+a[7]*x**7+a[8]*x**8 for x in np.arange(0,0.07,0.001)])
plt.xlabel('radial coordinate, m')
plt.ylabel('Intensity, W/m²')
plt.show()

a = a[1:]#we don't take a0 since we take the derivative of the polynome
myemeas = np.matmul(AbelMatrix,a)/-3.141592
#plt.plot(r,emeas)
#plt.show()

abel_obj = abel1D.AbelInversion1D(profile)
profile_inverted = abel_obj.invert()
emeas = profile_inverted.intensity
r = profile_inverted.r.value
while len(myemeas) < len(emeas):
    myemeas = np.append(myemeas,0)
plt.plot(r,emeas)
plt.plot(r,myemeas)
plt.xlabel('radial coordinate, mm')
plt.ylabel('Emission, W/m³')
plt.show()
#emeas=np.append(emeas,emeas)
#r=np.append(r,-r)
#plt.plot(r,emeas)
#plt.show()
end = 500
a = np.polyfit(r[0:end],emeas[0:end].value,7)
end = 500
#plt.plot(r[0:end],[a[7]+a[6]*x+a[5]*x**2+a[4]*x**3+a[3]*x**4+a[2]*x**5+a[1]*x**6+a[0]*x**7 for x in r[0:end]])
#plt.plot([r[end],90],[a[7]+a[6]*r[end]+a[5]*r[end]**2+a[4]*r[end]**3+a[3]*r[end]**4+a[2]*r[end]**5+a[1]*r[end]**6+a[0]*r[end]**7,0])
#plt.show()

mypts = [0,2.5,5,7.5,10,15,20,25,30,40]
#mypts = [0]
x = []
for p in mypts:
    x.append(min(range(len(r)), key=lambda i: abs(r[i]-p)))
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
    #nframe = int(random.uniform(1,11))
    #ncalibframe = int(random.uniform(1,6))
    #nframe = int(i/(nTest/10))+1 #take an other image every 1/10 of ntest
    nframe = -1# it means we generate a random image from mean and std deviation
    ncalibframe = 0
    L700_910 = meas700_910.calibration('randcorr', nframe, ncalibframe, mu, sigma)
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
plt.ylabel('temperature, K')
plt.show()

plt.hist(radialT[0],bins=20,density=True)
plt.xlabel('temperature, K')
plt.ylabel('Probability distribution, []')
plt.show()

plt.hist(radialT[9],bins=20,density=True)
plt.xlabel('temperature, K')
plt.ylabel('Probability distribution, []')
plt.show()
