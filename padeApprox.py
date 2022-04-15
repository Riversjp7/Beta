#!/usr/bin/env python3

import numpy as np
import pandas as pd
'''Pade Approximant'''

'''The coefficients of a are known and are constrained by L and M'''
L=2
M=2
aMax=L+M

'''Creation of an array of coefficients for e^x (Taylor series)'''
aPop=np.ones(aMax)
aMed=np.add.accumulate(aPop)
aMed=np.multiply.accumulate(aMed)
aMed=np.insert(aMed,0,1)
a=1/aMed
#print(aMed)

'''populate a coefficient array for P/Q based on L and M. 
   This creates array columns q1,...,qm,p1,...,pl'''
cofA=pd.DataFrame(np.zeros([aMax,aMax]))
for i in range(M):
    cofA.loc[i:, [i]] = a[:aMax-i,np.newaxis]
for j in range(L):
    cofA.loc[j, [M+j]] = -1
#print(cofPop,cofPop.loc[1:, [1]],a[:aMax-1,np.newaxis],sep='\n')
'''Now that that madness is done, time to solve the system of equations
   cofA*pqVal=-a, where cofA is the coefficients and X is the 
   array [q1,q2,...,qM,p1,p2,...pL]. 
   Note: a starts at index 1 not 0.'''

pqVal=np.linalg.solve(cofA,-a[1:])
#print(pqVal)

'''Now initial values and P/Q polynomials'''
xval=2
q0=1
p0=a[0]
Q=pqVal[:M]
P=pqVal[M:]
def pol(a0,pn,x):
    pn=np.insert(pn,0,a0)
    print(pn)
    xInit=np.full(len(pn)-1,x)
    xPol=np.multiply.accumulate(xInit)
    return a0+np.sum(pn[1:]*xPol)
pVal=pol(p0,P,xval)
qVal=pol(q0,Q,xval)


#print(pqVal, pVal/qVal,sep='\n')