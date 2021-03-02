from function import efromT
import random

#AAA ≤ 0.3%
#AA ≤ 1%
#A+ ≤ 2%
#A  ≤ 3%
#B+ ≤ 7%
#B  ≤ 10%
#C+ ≤ 18%
#C  ≤ 25%
#D+ ≤ 40%
#D  ≤ 50%
mtoJ = 1.986445857E-25

T = 7000 #[K]
P = 10000
Eu = [8663145.4*mtoJ,8662777.8*mtoJ,8662575.7*mtoJ] #from [m^-1] to [J]
El = [7376820*mtoJ,7376820*mtoJ,7376820*mtoJ]
g = [7,5,3]
Au = [3.69E7,3.69E7,3.69E7]

'''err on P'''
# +- 10% with constant probability distribution
maxerr = 10 #in %
Pmin = P - P/100* maxerr
Pmax = P + P/100* maxerr
Pvect = []
resolution = 1000
for i in range(resolution+1):
    Pvect.append(Pmin + i*(Pmax-Pmin)/resolution)
probdistr = []
for i in range(resolution+1):
    probdistr.append(1) #constant ditribution
total = sum(probdistr)

numtest = 1000
results = []
for j in range(numtest):
    x = random.uniform(0, total)
    i = -1
    while x > 0:
        i = i + 1
        x = x - probdistr[i]
    P = Pvect[i]
    e = efromT(T, P, 'O', g, Eu, Au, El)
    results.append(e)
    #print(j)

#print(results)
rangeresults = []
resolution = 10
for i in range(resolution+1):
    rangeresults.append(min(results) + i*(max(results)-min(results))/resolution)
#print(rangeresults)
resultsdistr = [0 for i in range(resolution)]
for e in results:
    i = 1
    while e > rangeresults[i]:
        i = i + 1
    resultsdistr[i-1] = resultsdistr[i-1] + 1/numtest
print(resultsdistr)
err = (max(results)-min(results))/(max(results)+min(results)*100)#in %
print(err)
