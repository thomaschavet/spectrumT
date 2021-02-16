import csv
import mpmath
import math

wl = []
I = []
with open('ps100mbar_T7000K_lambda300-900nm.csv') as csvfile:
    data = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    for row in data:
        wl.append(row[0])
        I.append(row[1])
wl.pop(0)
I.pop(0)

'''Oxygen I'''
start = 776.5 #[nm]
end = 778.0 #[nm]

#find closest to start and end
i = 0
while wl[i] < start:
    i = i + 1
#wl(i) >= start
if wl[i] - start < start - wl[i-1]:
    start = wl[i]
    starti = i
else:
    start = wl[i-1]
    starti = i-1
while wl[i] < end:
    i = i + 1
#wl(i) >= end
if wl[i] - end < end - wl[i-1]:
    end = wl[i]
    endi = i
else:
    end = wl[i-1]
    endi = i-1

emeas = 0
for i in range(starti, endi):
    baseline1 = I[starti] + (wl[i]-start)/(end-start)*(I[endi]-I[starti])
    baseline2 = I[starti] + (wl[i+1]-start)/(end-start)*(I[endi]-I[starti])
    emeas = emeas + ((I[i]-baseline1 + I[i+1]-baseline2)*(wl[i+1]-wl[i])/2) #(B+b)*h/2


P = 10000 #[pa]
Eu1 = 8663145.4 * 1.986445857E-25 #from [m^-1] to [J]
Eu2 = 8662777.8 * 1.986445857E-25
Eu3 = 8662575.7 * 1.986445857E-25
El = 7376820 * 1.986445857E-25
g1 = 7
g2 = 5
g3 = 3
gg = 5
Kb = 1.38064852E-23
Au1 = 3.69E7
Au2 = 3.69E7
Au3 = 3.69E7
#Guess
T = 7000 #[K]

err = 1000
while abs(err) > 1E-10:
    ng = P/(Kb*T)*0.2
    Qint = 5 + 3*mpmath.exp(-15826.5*1.986445857E-25/(Kb*T)) + 1*mpmath.exp(-22697.7*1.986445857E-25/(Kb*T)) + 5*mpmath.exp(-1586786.2*1.986445857E-25/(Kb*T)) + 1*mpmath.exp(-3379258.3*1.986445857E-25/(Kb*T)) + 5*mpmath.exp(-7376820.0*1.986445857E-25/(Kb*T)) 
    #print(Qint)
    nu1 = ng*g1*mpmath.exp(-Eu1/(Kb*T))/Qint
    nu2 = ng*g2*mpmath.exp(-Eu2/(Kb*T))/Qint
    nu3 = ng*g3*mpmath.exp(-Eu3/(Kb*T))/Qint
    e1 = (Eu1-El)/(4*math.pi) * Au1*nu1
    e2 = (Eu2-El)/(4*math.pi) * Au2*nu2
    e3 = (Eu3-El)/(4*math.pi) * Au3*nu3
    ecalc = e1 + e2 + e3
    err = emeas - ecalc
    eta = 0.1
    T = T + eta*err

#print(emeas)
#print(ecalc)
print('Oxygen I:')
print(T)
