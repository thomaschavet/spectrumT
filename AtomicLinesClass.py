#defines the transition of the atom. inputs are: type of atom ('O' or 'N'), energy of high level,
#energy of low level, Einstein coefficien and degeneracy level g.
#Output is total emission of that atom transition from temperature and pressure
import mpmath
import math
import sys

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
        #generate Au from constant probability distribution
        if type(self.AuVar) == type([]):
            self.Au = [self.AuVar[i].rand() for i in range(len(self.AuVar))]
        else:
            myAu= self.AuVar.rand()
            self.Au = [myAu for i in range(len(self.g))]
    
    def randAuGauss(self):
        #generate Au from Gaussian probability distribution
        if type(self.AuVar) == type([]):
            self.Au = [self.AuVar[i].randgauss() for i in range(len(self.AuVar))]
        else:
            myAu= self.AuVar.randgauss()
            self.Au = [myAu for i in range(len(self.g))]
    
    def emission(self, T, P):
        #search for the population value in 2D table computed with mutationpp
        file = open("PopTab.txt","r")
        nTab = file.read().splitlines()
        file.close()
        Pup = float(nTab[0].split()[1])
        i = 1
        while Pup < P:
            i = i + 1
            Pup = float(nTab[0].split()[i])
        Pdown = float(nTab[0].split()[i-1])
        ratioP = (P-Pup)/(Pdown-Pup)
        Tup = float(nTab[1].split()[0])
        j = 1
        while Tup < T:
            j = j + 1
            Tup = float(nTab[j].split()[0])
        Tdown = float(nTab[j-1].split()[0])
        ratioT = (T-Tup)/(Tdown-Tup)
        nup = float(nTab[j].split()[i])*(1-ratioT) + float(nTab[j-1].split()[i])*ratioT
        ndown = float(nTab[j].split()[i-1])*(1-ratioT) + float(nTab[j-1].split()[i-1])*ratioT
        n = nup*(1-ratioP) + ndown*ratioP
        
        if self.atom == 'O':
            file = open("PopOTab.txt","r")
            OTab = file.read().splitlines()
            file.close()
            Oup = float(OTab[j].split()[i])*(1-ratioT) + float(OTab[j-1].split()[i])*ratioT
            Odown = float(OTab[j].split()[i-1])*(1-ratioT) + float(OTab[j-1].split()[i-1])*ratioT
            XO = Oup*(1-ratioP) + Odown*ratioP
            ng = n*XO
        elif self.atom == 'N':
            file = open("PopNTab.txt","r")
            NTab = file.read().splitlines()
            file.close()
            Nup = float(NTab[j].split()[i])*(1-ratioT) + float(NTab[j-1].split()[i])*ratioT
            Ndown = float(NTab[j].split()[i-1])*(1-ratioT) + float(NTab[j-1].split()[i-1])*ratioT
            XN = Nup*(1-ratioP) + Ndown*ratioP
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
        
        #computes emission from the two equations
        Kb = 1.38064852E-23
        ecalc = 0
        for i in range(len(self.g)):
            nu = ng*self.g[i]*mpmath.exp(-self.Eu[i]/(Kb*T))/Qint
            e = (self.Eu[i]-self.El[i])/(4*math.pi) * self.Au[i]*nu
            ecalc = ecalc + e
        return ecalc