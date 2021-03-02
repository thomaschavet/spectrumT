import mpmath
import math
import os
import sys
import pathlib

def efromT(T, P, atom_spec, g, Eu, Au, El):
    path = pathlib.Path(__file__).parent.absolute()
    numfract = os.popen('''export MPP_DIRECTORY='''+str(path)+'''/Mutationpp
          export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
          export PATH=$MPP_DIRECTORY/install/bin:$PATH
          export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH
          '''+str(path)+'''/Mutationpp/src/general/mppequil --no-header -T ''' + str(T) + ''' -P ''' + str(P) + ''' -m 4 -s 0 plasmatron''').read().split()
    n = float(numfract[0])
    if atom_spec == 'O':
        XO = float(numfract[12])
        ng = n*XO
        partfunc = open("O.res","r")
        Qtab = partfunc.read().splitlines()
        partfunc.close()
    elif atom_spec == 'N':
        XN = float(numfract[11])
        ng = n*XN
        partfunc = open("N.res","r")
        Qtab = partfunc.read().splitlines()
        partfunc.close()
    else:
        sys.exit("atom_spec must be 'O' or 'N' !")
    Tdown = float(Qtab[0].split()[0])
    i = 0
    while Tdown > T:
        i = i + 1
        Tdown = float(Qtab[i].split()[0])
    Tup = float(Qtab[i-1].split()[0])
    ratio = (T-Tdown)/(Tup-Tdown)
    Qint = float(Qtab[i].split()[1])*(1-ratio) + float(Qtab[i-1].split()[1])*ratio
    Kb = 1.38064852E-23
    ecalc = 0
    for i in range(len(g)):
        nu = ng*g[i]*mpmath.exp(-Eu[i]/(Kb*T))/Qint
        e = (Eu[i]-El[i])/(4*math.pi) * Au[i]*nu
        ecalc = ecalc + e
    return ecalc