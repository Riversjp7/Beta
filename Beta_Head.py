#!/usr/bin/env python3
import numpy as np
import pandas as pd
import numpy.polynomial.polynomial as poly
from matplotlib import pyplot as plt
import tkinter as tk
from matplotlib.widgets import Button, Slider

'''This is a beta object complete with Pade approx. capability'''
class Beta:
	def __init__(self, nf=12,L=3,M=3,a=0.634/np.pi,a_range=None):
		self.nf=nf
		self.M=M
		self.L=L
		self.a=a
		self.aMax=L+M
		if a_range is None:
			a_range=np.linspace(0,1,100)
		self.a_range=a_range
		zeta=np.array([-0.5,0,
				   1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587053040043189623371906796287246,
				   1.2020569031595942853997381615114499907649862923404988817922715553,
				   1.0823232337111381915160036965411679027747509519187269076829762154,
				   1.0369277551433699263313654864570341680570809195019128119741926779])
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
				 
		self.bCo=-1*np.array([0,0,b0,b1,b2,b3,b4])
	def set_nf(self, nf):
		self.nf=nf
	def set_LM(self, L,M):
		self.L=L
		self.M=M
		self.aMax=M+L
	def set_a_range(self, a_Max,a_Min,a_Step):
		self.a_range=np.linspace(a_Max,a_Min,a_Step)
	def set_a(self,a):
		self.a=a
	def pade_solve(self):
		'''This function use a matrix to solve the system of eqs. It returns and array
		of coefficients for p(x) and q(x)'''
		cofA=pd.DataFrame(np.zeros([self.aMax,self.aMax]))
		for i in range(self.M):
			cofA.loc[i:, [i]] = self.bCo[:self.aMax-i,np.newaxis]
		for j in range(self.L):
			cofA.loc[j, [self.M+j]] = -1
		pqVal=np.linalg.solve(cofA,self.bCo[1:])
		q_cof=pqVal[:self.M]
		p_cof=pqVal[self.M:]
		q_cof=np.insert(q_cof,0,1)
		p_cof=np.insert(p_cof,0,self.bCo[0])
		return p_cof,q_cof 
	def pade_val(self):
		'''This function uses the coefficient arrays to compute an actual value for
		p/q'''
		pcof,qcof=self.pade_solve()
		xpInit=np.full(len(pcof)-1,self.a)
		xqInit=np.full(len(qcof)-1,self.a)
		xPPol=np.multiply.accumulate(xpInit)
		xqPol=np.multiply.accumulate(xqInit)
		p_over_q=(pcof[0]+np.sum(pcof[1:]*xPPol))/(qcof[0]+np.sum(qcof[1:]*xqPol))
		return p_over_q
	def curve(self,pn=None,x=None):
		'''Here we take polynomial coefficients and add a range of x values for corresponding orders
		constant term is added and the input coefficient array is mapped to increasing orders of x'''
		if pn is None:
			pn=self.bCo
		self.pn=pn
		if x is None:
			x=self.a_range[:,np.newaxis]
		xInit=np.full((len(x),len(pn)-1),x)
		xPol=np.multiply.accumulate(xInit,1)
		padPol=pn[0]+np.sum(pn[1:]*xPol,axis=1)
		return padPol
	def zeros(self,pn=None):
		if pn is None:
			pn=self.curve()
		self.pn=pn
		return poly.polyroots(self.pn)