import csv


def efromSpectr(start, end):
    wl = []
    I = []
    with open('ps100mbar_T7000K_lambda300-900nm.csv') as csvfile:
        data = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in data:
            wl.append(row[0])
            I.append(row[1])
    wl.pop(0)
    I.pop(0)
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


def Tcalc(e, atomLines, P):
    #Guess
    T = 7000 #[K]
    err = 1000
    newerr = 0
    eta = 0.8
    while abs(err) > 1E-3:
        ecalc = atomLines.emission(T, P)
        newerr = e - ecalc
        if (newerr>=0 and err>=0) or (newerr<0 and err<0):
            eta = eta + eta*0.5
        else:
            eta = eta - eta*0.2
        err = newerr
        T = T + eta*err
    return T

def fact(n):
    fact = 1
    for num in range(2, n + 1):
        fact *= num
    return fact

def integral(x,r,m):#solve integral of m*x^(m-1)/(sqrt(x²-r²))
    #return m*x**(m-1)/(x**2-r**2)**(1/2)
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
            #newline.append(scipy.integrate.quad(integral,r[i],R,args=(r[i],m))[0])
        matrix.append(newline)
    return matrix
