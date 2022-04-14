#!/usr/bin/env python3

from Beta_Head import *

b=Beta()

#some initial conditions for the plots
scaleFac=12
fig, ax =plt.subplots(1,3,figsize=[scaleFac,scaleFac*9/16])
fig.tight_layout(pad=3.0)
plt.subplots_adjust(bottom=0.25,wspace=0.25)
	
def asymp(pade):
	'''This is a small function to eliminate Pade asymptotes connected continuously in 
	the graph. It doesn't work well.''''
	ll=-0.1
	hl=0.1
	pade[pade>hl]=np.inf
	pade[pade<ll]=-np.inf
	return pade
def on_graph(b):
	'''This function is called upon by the slider and graphs the Pade approximates 
	appropriately'''
	for i in range(3):
		ax[i].clear()
		ax[i].plot(b.a_range, 0*b.a_range,'k--')
		ax[i].set_ylim([-0.05,0.05])
		ax[i].set_xlabel(r'coupling constant $a_s$')
		ax[i].set_ylabel(r'$\beta(a)$')
		ax[i].set_title(r'$\beta(a)$ loop 5 and Pade Approx.')
		ax[i].plot(b.a_range,b.curve(),label=f'loop order 5')
	p,q=b.pade_solve()
	ax[2].plot(b.a_range,asymp(-b.curve(p)/b.curve(q)),label=r'$\beta_{3,3}$')
	b.set_LM(4,2)
	p,q=b.pade_solve()
	ax[1].plot(b.a_range,asymp(-b.curve(p)/b.curve(q)),label=r'$\beta_{4,2}$')
	b.set_LM(2,4)
	p,q=b.pade_solve()
	ax[0].plot(b.a_range,asymp(-b.curve(p)/b.curve(q)),label=r'$\beta_{2,4}$')
	for i in range(3):
		ax[i].legend()
	
on_graph(b) # initial graph

'''Slider functionality'''
nfAx=plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='lightgrey')
nfSlide = Slider(nfAx, 'Flavor # ', 0, 17, valinit=b.nf, valstep=1)
	
def update_nf(val):
	flav = nfSlide.val
	beta=Beta(nf=flav)
	on_graph(beta)
	fig.canvas.draw_idle()
nfSlide.on_changed(update_nf)
plt.show()
