#!/usr/bin/env python3

'''Pade Approximant of beta function'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox



#array of zeta values (1 is 0 and not inf)
zeta=np.array([-0.5,0,
               1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587053040043189623371906796287246,
               1.2020569031595942853997381615114499907649862923404988817922715553,
               1.0823232337111381915160036965411679027747509519187269076829762154,
               1.0369277551433699263313654864570341680570809195019128119741926779])
nf=12
b0=(0.25)*(11-2*nf/3)
b1=(0.25**2)*(102-38*nf/3)
b2=(0.25**3)*(2857/2-5033*nf/18+325*(nf**2)/54)
b3=(0.25**4)*(149753/6+3564*zeta[3]-(1078361/162+6508*zeta[3]/27)*nf
				 +(50056/162+6472*zeta[3]/81)*nf**2
				 +(1093/729)*nf**3)
b4=(0.25**5)*(8157455/16+621885*zeta[3]/2-88209*zeta[4]/2-288090*zeta[5]
				 +(-336460813/1944-4811164*zeta[3]/81+33935*zeta[4]/6 + 1358995*zeta[5]/27)*nf
				 +(25960913/1944+698531*zeta[3]/81-10526*zeta[4]/9-381760*zeta[5]/81)*nf**2
				 +(-630559/5832-48722*zeta[3]/243+1618*zeta[4]/27+460*zeta[5]/9)*nf**3
				 +(1205/2916-152*zeta[3]/81)*nf**4)
				 
bCo=np.array([0,0,b0,b1,b2,b3,b4])
'''populate a coefficient array for P/Q based on L and M. 
   This creates array columns q1,...,qm,p1,...,pl'''
'''The coefficients of a are known and are constrained by L and M'''
M=3
L=3
aMax=M+L
cofA=pd.DataFrame(np.zeros([aMax,aMax]))
for i in range(M):
    cofA.loc[i:, [i]] = bCo[:aMax-i,np.newaxis]
for j in range(L):
    cofA.loc[j, [M+j]] = -1
    
pqVal=np.linalg.solve(cofA,bCo[1:])


'''Now initial values and P/Q polynomials'''
xtest=2
xval=[0.129711278619895,0.634/np.pi,0.614/np.pi]
q0=1
p0=bCo[0]
Q=pqVal[:M]
P=pqVal[M:]
def pol(a0,pn,x):
    #pn=np.insert(pn,0,a0)
    #print(pn,x,len(pn),sep='\n')
    xInit=np.full((len(x),len(pn)),x)
    xPol=np.multiply.accumulate(xInit,0)
    padPol=a0+np.sum(pn[:,np.newaxis]*xPol,axis=0)
    print()
    #print(padPol,pn[:,np.newaxis]*xPol,sep='\n')
    return padPol
pVal=pol(p0,P,xtest)
#qVal=pol(q0,Q,xval)

#print(p0,pVal,'\n',q0,qVal)
#print(pVal/qVal,sep='\n')