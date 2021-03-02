import math
from function import efromT

#CONSTANT:
mtoJ = 1.986445857E-25
Kb = 1.38064852E-23
P = 10000

'''Oxygen I'''
Eu = [8663145.4*mtoJ,8662777.8*mtoJ,8662575.7*mtoJ] #from [m^-1] to [J]
El = [7376820*mtoJ,7376820*mtoJ,7376820*mtoJ]
g = [7,5,3]
Au = [3.69E7,3.69E7,3.69E7]

meanT = 7000
allA = [0,1,2,3,4,5,6,7,8,9,10] #Amplitude of the oscilation of T in %
error = []
for A in allA:
    DeltT = A/100*meanT
    step = 20
    emean = 0
    for i in range(step):
        T = meanT + math.sin(i/step*2*math.pi)*DeltT
        emean = emean + efromT(T, P, 'O', g, Eu, Au, El)
    emean = emean/step
    etrue = efromT(meanT, P, 'O', g, Eu, Au, El)
    #guess
    T = 7000
    err = 1000
    while abs(err) > 1E-5:
        ecalc = efromT(T, P, 'O', g, Eu, Au, El)
        err = emean - ecalc
        eta = 0.5
        T = T + eta*err
    print("err of T (in %) with A = "+ str(A) +":")
    print((T-meanT)/meanT*100)
    error.append(str((T-meanT)/meanT*100))
output = open("outputOI.txt","w")
for i in range(len(error)):
    output.write(error[i]+'\n')
output.close()

'''Oxygen II'''
Eu = [8863130.3*mtoJ,8863114.6*mtoJ,8863058.7*mtoJ] #from [m^-1] to [J]
El = [7679497.8*mtoJ,7679497.8*mtoJ,7679497.8*mtoJ]
g = [1,5,3]
Au = [3.22E7,3.22E7,3.22E7]

meanT = 7000
allA = [0,1,2,3,4,5,6,7,8,9,10] #Amplitude of the oscilation of T in %
error = []
for A in allA:
    DeltT = A/100*meanT
    step = 20
    emean = 0
    for i in range(step):
        T = meanT + math.sin(i/step*2*math.pi)*DeltT
        emean = emean + efromT(T, P, 'O', g, Eu, Au, El)
    emean = emean/step
    etrue = efromT(meanT, P, 'O', g, Eu, Au, El)
    #guess
    T = 7000
    err = 1000
    while abs(err) > 1E-5:
        ecalc = efromT(T, P, 'O', g, Eu, Au, El)
        err = emean - ecalc
        eta = 0.5
        T = T + eta*err
    print("err of T (in %) with A = "+ str(A) +":")
    print((T-meanT)/meanT*100)
    error.append(str((T-meanT)/meanT*100))
output = open("outputOII.txt","w")
for i in range(len(error)):
    output.write(error[i]+'\n')
output.close()

'''Nitrogen I'''
Eu = [9675084 * mtoJ] #from [m^-1] to [J]
El = [8336462 * mtoJ]
g = [4]
Au = [1.96E7]

meanT = 7000
allA = [0,1,2,3,4,5,6,7,8,9,10] #Amplitude of the oscilation of T in %
error = []
for A in allA:
    DeltT = A/100*meanT
    step = 20
    emean = 0
    for i in range(step):
        T = meanT + math.sin(i/step*2*math.pi)*DeltT
        emean = emean + efromT(T, P, 'N', g, Eu, Au, El)
    emean = emean/step
    etrue = efromT(meanT, P, 'N', g, Eu, Au, El)
    #guess
    T = 7000
    err = 1000
    while abs(err) > 1E-5:
        ecalc = efromT(T, P, 'N', g, Eu, Au, El)
        err = emean - ecalc
        eta = 0.5
        T = T + eta*err
    print("err of T (in %) with A = "+ str(A) +":")
    print((T-meanT)/meanT*100)
    error.append(str((T-meanT)/meanT*100))
output = open("outputNI.txt","w")
for i in range(len(error)):
    output.write(error[i]+'\n')
output.close()

'''Nitrogen II'''
Eu = [9488182*mtoJ,9483089*mtoJ,9479349*mtoJ] #from [m^-1] to [J]
El = [8336462*mtoJ,8331783*mtoJ,8328407*mtoJ]
g = [8,6,4]
Au = [2.53E7,1.88E7,1.15E7]

meanT = 7000
allA = [0,1,2,3,4,5,6,7,8,9,10] #Amplitude of the oscilation of T in %
error = []
for A in allA:
    DeltT = A/100*meanT
    step = 20
    emean = 0
    for i in range(step):
        T = meanT + math.sin(i/step*2*math.pi)*DeltT
        emean = emean + efromT(T, P, 'N', g, Eu, Au, El)
    emean = emean/step
    etrue = efromT(meanT, P, 'N', g, Eu, Au, El)
    #guess
    T = 7000
    err = 1000
    while abs(err) > 1E-3:
        ecalc = efromT(T, P, 'N', g, Eu, Au, El)
        err = emean - ecalc
        eta = 0.5
        T = T + eta*err
    print("err of T (in %) with A = "+ str(A) +":")
    print((T-meanT)/meanT*100)
    error.append(str((T-meanT)/meanT*100))
output = open("outputNII.txt","w")
for i in range(len(error)):
    output.write(error[i]+'\n')
output.close()
