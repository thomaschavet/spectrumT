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
    while abs(err) > 1E-3:
        ecalc = atomLines.emission(T, P)
        err = e - ecalc
        eta = 0.8
        T = T + eta*err
    return T
