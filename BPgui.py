#!/usr/bin/env python3

'''Pade Approximant of beta function'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import tkinter as tk
import numpy.polynomial.polynomial as poly
import subprocess as sp
#import os

def pol(a0,pn,x):
    pn=np.insert(pn,0,a0)
    #print(pn*120)
    xInit=np.full(len(pn)-1,x)
    xPol=np.multiply.accumulate(xInit)
    return a0+np.sum(pn[1:]*xPol),pn
def on_push():
	#array of zeta values (1 is 0 and not inf)
	zeta=np.array([-0.5,0,
				   1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587053040043189623371906796287246,
				   1.2020569031595942853997381615114499907649862923404988817922715553,
				   1.0823232337111381915160036965411679027747509519187269076829762154,
				   1.0369277551433699263313654864570341680570809195019128119741926779])
	nf=float(ent_nf.get())
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
	M=int(ent_M.get())
	L=int(ent_L.get())
	aMax=M+L
	cofA=pd.DataFrame(np.zeros([aMax,aMax]))
	for i in range(M):
		cofA.loc[i:, [i]] = bCo[:aMax-i,np.newaxis]
	for j in range(L):
		cofA.loc[j, [M+j]] = -1
	
	pqVal=np.linalg.solve(cofA,bCo[1:])
	Q=pqVal[:M]
	P=pqVal[M:]
	print(Q,P,sep='\n')
	'''Now initial values and P/Q polynomials'''
	xval=float(ent_a.get())
	q0=1
	p0=bCo[0]


	pVal,pCof=pol(p0,P,xval)
	qVal,qCof=pol(q0,Q,xval)
	poverq=pVal/qVal
	lbl_L["text"] = f"L = {L}"
	lbl_M["text"] = f"M = {M}"
	lbl_a["text"] = f"a = {xval}"
	lbl_nf["text"]=f"nf = {nf}"
	lbl_pq["text"] = f"P/Q = {poverq}"
	lbl_p["text"] = f"Coeffiecients of p: {pCof[:]}"
	lbl_q["text"] = f"Coeffiecients of q:{qCof[:]}"
	lbl_zeros["text"] = f"Zeroes:{poly.polyroots(pCof)}"
	#print(poly.polyroots(pCof))
	
	return pVal,pCof,qVal,qCof
def graph_func():
	pv,pc,qv,qc=on_push()
	
	
	
	
	'''Creating a tkinter  win'''
entWidth=25
win=tk.Tk()
win.rowconfigure(0,minsize=25,weight=1)
win.columnconfigure([0,1],minsize=50,weight=1)
win.title(f'Pade Approx. for \u03b2 function')

frm_in=tk.Frame(master=win)
frm_in.grid(row=0,column=0)
frm_out=tk.Frame(master=win)
frm_out.grid(row=0,column=1)


ent_L=tk.Entry(master=frm_in,width=entWidth)
ent_L.grid(row=0,column=1,padx=10)
ent_L.insert(tk.END, 3)
ent_M=tk.Entry(master=frm_in,width=entWidth)
ent_M.grid(row=1,column=1,padx=10)
ent_M.insert(tk.END, 3)
ent_a=tk.Entry(master=frm_in,width=entWidth)
ent_a.grid(row=2,column=1,padx=10)
ent_a.insert(tk.END, 0.634/np.pi)
ent_nf=tk.Entry(master= frm_in,width=entWidth)
ent_nf.grid(row=3,column=1,padx=10)
ent_nf.insert(tk.END, 12)

lbl_L=tk.Label(master= frm_in,text=f'L value = {ent_L.get()} ')
lbl_L.grid(row=0,column=0,padx=10)
lbl_M=tk.Label(master= frm_in,text=f'M value = {ent_M.get()} ')
lbl_M.grid(row=1,column=0,padx=10)
lbl_a=tk.Label(master=frm_in,text=f'a value = {ent_a.get()} ')
lbl_a.grid(row=2,column=0,padx=10)
lbl_nf=tk.Label(master=frm_in,text=f'nf value = {ent_nf.get()} ')
lbl_nf.grid(row=3,column=0,padx=10)

#p,pc,q,qc=on_push()
btn_up=tk.Button(master= frm_in, text="Submit",command=on_push)
btn_up.grid(row=4,column=0,sticky="nsew")
btn_up=tk.Button(master= frm_in, text="Graph",command=on_push)
btn_up.grid(row=4,column=2,sticky="nsew")

lbl_pq=tk.Label(master= frm_out,text=f'P/Q value =  ')
lbl_pq.grid(row=0,column=2, sticky="w")
lbl_p=tk.Label(master= frm_out,text=f'Coefficients of P = ')
lbl_p.grid(row=1,column=2,sticky="w")
lbl_q=tk.Label(master= frm_out,text=f'Coefficients of Q=  ')
lbl_q.grid(row=2,column=2, sticky="w")

lbl_zeros=tk.Label(master= frm_out,text=f'zeros: ')
lbl_zeros.grid(row=3,column=2, sticky="w")
'''lbl_scale=tk.Label(master= frm_out,text=f'Scale value = none ')
lbl_scale.grid(row=3,column=2, sticky="e")

ent_scale=tk.Entry(master= frm_out,width=entWidth)
ent_scale.grid(row=3,column=3)

btn_scale=tk.Button(master= frm_out, text="Scale")
btn_scale.grid(row=4,column=3,sticky="nsew")'''

win.mainloop()