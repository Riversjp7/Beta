#!/usr/bin/env python3
'''This code plots the Beta function as a function of 
    quark-gluon coupling with the number of flavors as a free parameter'''

import numpy as np
from matplotlib.widgets import Button, Slider
import matplotlib.pyplot as plt

fig, ax =plt.subplots()
plt.subplots_adjust(left=0.25,bottom=0.25)
#number of flavors
nf=6
nf0=6
a=np.arange(0,1,0.0001)
alph=a[:,np.newaxis]
p=[]


#array of zeta values (1 is 0 and not inf)
zeta=np.array([-0.5,0,
               1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587053040043189623371906796287246,
               1.2020569031595942853997381615114499907649862923404988817922715553,
               1.0823232337111381915160036965411679027747509519187269076829762154,
               1.0369277551433699263313654864570341680570809195019128119741926779])

# Two of the beta functions that include zeta functions
b0=(1/4)*(11-2*nf/3)
b1=(1/4**2)*(102-18*nf/3)
b2=(1/4**3)*(2857/2-5033*nf/18+325*(nf**2)/54)
b3=(1/4**4)*(149753/6+3564*zeta[3]-(1078361/162+6508*zeta[3]/27)*nf
             +(50056/162+6472*zeta[3]/81)*nf**2
             +(1093/729)*nf**3)
b4=(1/4**5)*(8157455/16+621885*zeta[3]/2-88209*zeta[4]/2-288090*zeta[5]
             +(-336460813/1944-4811164*zeta[3]/81+33935*zeta[4]/6 + 1358995*zeta[5]/27)*nf
             +(25960913/1944+698531*zeta[3]/81-10526*zeta[4]/9-381760*zeta[5]/81)*nf**2
             +(-630559/5832-48722*zeta[3]/243+1618*zeta[4]/27+460*zeta[5]/9)*nf**3
             +(1205/2916-152*zeta[3]/81)*nf**4)
#pre-reduced functions
betRed3=114.23-27.1339*nf+1.58238*nf**2+0.0058567*nf**3
betRed4=524.56-181.8*nf+17.16*nf**2-0.22586*nf**3-0.0017933*nf**4

'''one array has coefficients, another array is made to have consequtive
polinomial terms (a^2, a^3...etc.), they are then multiplied and each
column of the resulting array is added to the index below it (this gives
the recursive relationship b_n=b_(n-1)+b_n*a^(n+2))'''
bCo=np.array([b0,b1,b2,b3,b4])
alphRay=np.full((len(alph),len(bCo)),alph)
alphPow=np.multiply.accumulate(alphRay, 1)
alphPow=alphPow*alph
bTrueCo=alphPow*bCo
bAcc=np.add.accumulate(bTrueCo,1)
plt.xlabel('coupling a')
plt.ylabel(r'$\beta$ term')
for i in range(len(bAcc[0,:])):
    p.append(plt.plot(alph,-bAcc[:,i],label=f'beta {i}'))
plt.legend()
plt.ylim([-10,10])

#nf slider
nfAx=plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='lightgrey')
nfSlide = Slider(nfAx, 'Flavor # ', 0, 30.0, valinit=nf0, valstep=1)
def update_nf(val):
	global nf
	nf = nfSlide.val
	fig.canvas.draw_idle()
# l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    #fig.canvas.draw_idle()

nfSlide.on_changed(update_nf)
print(nf,sep='\n')
plt.show()

