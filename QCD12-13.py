#!/usr/bin/env python3
'''This code plots the Beta function as a function of 
    quark-gluon coupling with the number of flavors as a free parameter'''

import numpy as np
from matplotlib.widgets import Button, Slider
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly

#some initial conditions for the plots
scaleFac=12
fig, ax =plt.subplots(1,3,figsize=[scaleFac,scaleFac*9/16])

fig.tight_layout(pad=3.0)
plt.subplots_adjust(bottom=0.25,wspace=0.25)

#number of flavors
nf0=12.8
a=np.arange(0,1,0.01)
alph0=a[:,np.newaxis]


def draw(y,Roots):
	'''This function is used to draw three graphs: B vs a, complex roots, list of complex roots.'''
	for i in range(3):
		ax[i].clear()
	ax[2].set_ylim([-0.05,0.05])
	ax[2].set_xlabel(r'coupling constant $a_s$')
	ax[2].set_ylabel(r'$\beta(a)$')
	ax[2].set_title(r'$\beta(a)$ function to loop order 5')
	for i in range(len(y[0,:])):
		ax[2].plot(alph0,-y[:,i],label=f'beta {i}')
	plt.legend()

	#plotting zeros
	ax[1].set_xlabel(r'real')
	ax[1].set_ylabel(r'imaginary')
	ax[1].set_title(r'$\beta$ roots in the complex plane')
	ax[1].scatter(Roots.real,Roots.imag,marker='o')

	#plotting telemetry
	ax[0].axis('off')
	ax[0].set_title('Complex Root Telemetry')
	bRootc=Roots[:,np.newaxis]
	for i in range(len(bRoot)):
		PForm=f'{np.around(bRootc[i].real,decimals=3)} + {np.around(bRootc[i].imag,decimals=3)}j'
		ax[0].text(0.5,0.9-0.1*i,PForm)
		
#array of zeta values (1 is 0 and not inf)
zeta=np.array([-0.5,0,
               1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587053040043189623371906796287246,
               1.2020569031595942853997381615114499907649862923404988817922715553,
               1.0823232337111381915160036965411679027747509519187269076829762154,
               1.0369277551433699263313654864570341680570809195019128119741926779])

# Two of the beta functions that include zeta functions
def beta(alph, nf):
	'''This function takes in a linspace array and a value for flavors and returns beta function values
	the roots of the beta function'''
	#coefficients of beta
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
	'''Here I am inserting 0s in a different array for finding roots'''
	bAt0=np.insert(bCo,[0,0],0)
	
	return bAcc,poly.polyroots(bAt0)

#plot beta vs a
yval,bRoot=beta(alph0,nf0)
draw(yval,bRoot)


#nf slider
nfAx=plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='lightgrey')
nfSlide = Slider(nfAx, 'Flavor # ', 12, 13, valinit=nf0, valstep=0.001)
def update_nf(val):
	flav = nfSlide.val
	newY,newRoot=beta(alph0, flav)
	draw(newY,newRoot)
	fig.canvas.draw_idle()
nfSlide.on_changed(update_nf)


plt.show()

