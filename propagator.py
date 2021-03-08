from RandVarClass import randVar
from AtomicLinesClass import atomLines
from function import efromSpectr
from function import Tcalc
import matplotlib.pyplot as plt

#Definition of the atomic lines (Here OI between 776 and 778 nm)
mtoJ = 1.986445857E-25
Eu = [8663145.4*mtoJ,8662777.8*mtoJ,8662575.7*mtoJ] #from [m^-1] to [J]
El = 7376820*mtoJ
g = [7,5,3]
meanAu = 3.69E7
Au = randVar(meanAu,3)
OIlines = atomLines('O', Eu, El, Au, g)

meanP = 10000 #[pa]
errP = 0.2 #error on the pressure in %
P = randVar(meanP, errP)

emeas = efromSpectr(776.5, 778.0)
e = randVar(emeas, 20)

results = []
nTest = 10000
for i in range(nTest):
    OIlines.randAu()#Take a new value for Au at each iteration
    T = Tcalc(e.rand(),OIlines,P.rand())
    results.append(float(T))
    print(i)

plt.hist(results,density=True)
plt.xlabel('Temperature')
plt.ylabel('probability distribution')
plt.show()