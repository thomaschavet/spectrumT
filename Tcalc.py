import csv
from function import efromT

wl = []
I = []
with open('ps100mbar_T7000K_lambda300-900nm.csv') as csvfile:
    data = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    for row in data:
        wl.append(row[0])
        I.append(row[1])
wl.pop(0)
I.pop(0)


def efromSpectr(start, end):
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
    return emeas

#CONSTANT:
mtoJ = 1.986445857E-25
Kb = 1.38064852E-23
P = 10000

'''Oxygen I'''

emeas = efromSpectr(776.5, 778.0)

Eu = [8663145.4*mtoJ,8662777.8*mtoJ,8662575.7*mtoJ] #from [m^-1] to [J]
El = [7376820*mtoJ,7376820*mtoJ,7376820*mtoJ]
g = [7,5,3]
Au = [3.69E7,3.69E7,3.69E7]
#Guess
T = 7000 #[K]

err = 1000
while abs(err) > 1E-5:
    ecalc = efromT(T, P, 'O', g, Eu, Au, El)
    err = emeas - ecalc
    eta = 1
    T = T + eta*err

print('Oxygen I:')
print(T)

'''Oxygen II'''
emeas = efromSpectr(844.2, 845.2)

Eu = [8863130.3*mtoJ,8863114.6*mtoJ,8863058.7*mtoJ] #from [m^-1] to [J]
El = [7679497.8*mtoJ,7679497.8*mtoJ,7679497.8*mtoJ]
g = [1,5,3]
Au = [3.22E7,3.22E7,3.22E7]

#Guess
T = 7000 #[K]

err = 1000
while abs(err) > 1E-5:
    ecalc = efromT(T, P, 'O', g, Eu, Au, El)
    err = emeas - ecalc
    eta = 1
    T = T + eta*err

print('Oxygen II:')
print(T)

'''Nitrogen I'''
emeas = efromSpectr(746.2, 747.4)

Eu = [9675084 * mtoJ] #from [m^-1] to [J]
El = [8336462 * mtoJ]
g = [4]
Au = [1.96E7]
#Guess
T = 7000 #[K]

err = 1000
while abs(err) > 1E-5:
    ecalc = efromT(T, P, 'N', g, Eu, Au, El)
    err = emeas - ecalc
    eta = 1
    T = T + eta*err

print('Nitrogen I:')
print(T)

'''Nitrogen II'''
emeas = efromSpectr(867.5, 869.1)

Eu = [9488182*mtoJ,9483089*mtoJ,9479349*mtoJ] #from [m^-1] to [J]
El = [8336462*mtoJ,8331783*mtoJ,8328407*mtoJ]
g = [8,6,4]
Au = [2.53E7,1.88E7,1.15E7]
#Guess
T = 7000 #[K]

err = 1000
while abs(err) > 1E-5:
    ecalc = efromT(T, P, 'N', g, Eu, Au, El)
    err = emeas - ecalc
    eta = 1
    T = T + eta*err

print('Nitrogen II:')
print(T)
