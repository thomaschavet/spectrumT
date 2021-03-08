import mpmath
import math
import os
import sys
import pathlib

class atomLines:
    #inputs must be [list] but can also be a single variable if it is the same for all lines
    def __init__(self, Atom, Eu, El, Au, g):
        self.atom = Atom
        self.AuVar = Au
        nlines = 1
        if type(Eu) == type([]):
            self.Eu = Eu
            nlines = len(Eu)
        if type(El) == type([]):
            self.El = El
            nlines = len(El)
        if type(g) == type([]):
            self.g = g
            nlines = len(g)
        
        #if input is not a list, we generate a list of that variable
        if type(Eu) != type([]):
            self.Eu = [Eu for i in range(nlines)]
        if type(El) != type([]):
            self.El = [El for i in range(nlines)]
        if type(g) != type([]):
            self.g = [g for i in range(nlines)]
        
        #we generate a random number for Au
        self.randAu()
        
        #store the tabulated value of Qint
        partfunc = open(Atom + ".res","r")
        self.Qtab = partfunc.read().splitlines()
        partfunc.close()
    
    def randAu(self):
        if type(self.AuVar) == type([]):
            self.Au = [self.AuVar[i].rand() for i in range(len(self.AuVar))]
        else:
            myAu= self.AuVar.rand()
            self.Au = [myAu for i in range(len(self.g))]
    
    def emission(self, T, P):
        #run mutationpp
        path = pathlib.Path(__file__).parent.absolute()
        numfract = os.popen('''export MPP_DIRECTORY='''+str(path)+'''/Mutationpp
              export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
              '''+str(path)+'''/Mutationpp/src/general/mppequil --no-header -T ''' + str(T) + ''' -P ''' + str(P) + ''' -m 4 -s 0 plasmatron''').read().split()
        n = float(numfract[0])
        if self.atom == 'O':
            XO = float(numfract[12])
            ng = n*XO
        elif self.atom == 'N':
            XN = float(numfract[11])
            ng = n*XN
        else:
            sys.exit("atom_spec must be 'O' or 'N' !")
        
        #computes Qint from tabulated values
        Tdown = float(self.Qtab[0].split()[0])
        i = 0
        while Tdown > T:
            i = i + 1
            Tdown = float(self.Qtab[i].split()[0])
        Tup = float(self.Qtab[i-1].split()[0])
        ratio = (T-Tdown)/(Tup-Tdown)
        Qint = float(self.Qtab[i].split()[1])*(1-ratio) + float(self.Qtab[i-1].split()[1])*ratio
        
        Kb = 1.38064852E-23
        ecalc = 0
        for i in range(len(self.g)):
            nu = ng*self.g[i]*mpmath.exp(-self.Eu[i]/(Kb*T))/Qint
            e = (self.Eu[i]-self.El[i])/(4*math.pi) * self.Au[i]*nu
            ecalc = ecalc + e
        return ecalc