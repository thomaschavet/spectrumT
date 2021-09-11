import math

def Tcalc(e, atomLines, P):
    #Guess
    T = 7000 #[K]
    err = 1000
    out = 0
    while abs(err) > 1E-3 and out == 0:
        #computes T with Newton-Raphson method
        ecalc = atomLines.emission(T, P)
        ecalc1 = atomLines.emission(T+1, P)
        slope = ecalc1 - ecalc
        err = e - ecalc
        T = T + err/slope
        if T < 2000:
            T=2000
            out = 1
    return T

def fact(n):
    fact = 1
    for num in range(2, n + 1):
        fact *= num
    return fact

def integral(x,r,m):#solve integral of m*x^(m-1)/(sqrt(x²-r²))
    if m%2==0:
        p = int(m/2-1)
        mysum = 0
        for k in range(p+1):
            mysum = mysum + fact(2*k)*fact(p)**2/(fact(2*p+1)*fact(k)**2)*(4*r**2)**(p-k)*x**(2*k)
        result = m*(x**2-r**2)**(1/2)*mysum
    else:
        p = int((m-1)/2)
        mysum = 0
        for k in range(1,p+1):
            mysum = mysum + fact(k)*fact(k-1)/fact(2*k)*r**(2*(p-k))*(2*x)**(2*k-1)
        result = m*fact(2*p)/(2**(2*p)*fact(p)**2)*((x**2-r**2)**(1/2)*mysum + r**(2*p)*math.log(x+(x**2-r**2)**(1/2)))
    return result
    
def AbelInvMatrix(r,deg):
    R = r[-1]
    matrix = []
    for i in range(len(r)):
        newline = []
        for m in range(1,deg+1):
            newline.append(integral(R,r[i],m)-integral(r[i],r[i],m))
        matrix.append(newline)
    return matrix