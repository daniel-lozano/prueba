import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad
from scipy.interpolate import interp1d

from sys import argv


FILE=open("t_points.dat","r")
data=FILE.readlines()


#Setting arrays for the calculus---------------------------------
f=[]
TP1=[]
TP2=[]
TP1_C=[]
TP2_C=[]


#Constants used-------------------------------------------------

gamma=float(data[1].split()[1])

alphaI=0.28125
alphaN=1.38
m=0
Io=0.904#
B=0.51/2

mu=1.0 #a.u.
hbar=1.0 #a.u.

h=2.0*np.pi
c=137.0
lam=735.0/0.0529

Z=2.0
omega=h*c/lam #
number=0.0#float(argv[1])

Factor=24.18884 #time[as]/a.u.


#Reading the important values----------------------------------
print("gamma=",gamma)

for i in range(1,len(data)):
    f.append(data[i].split()[0])
    TP1.append(data[i].split()[2])
    TP2.append(data[i].split()[3])
    TP1_C.append(data[i].split()[4])
    TP2_C.append(data[i].split()[5])
    


#------------------------------------------Energy factor-------------------------------------------
#D=0 #Dissipation

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)


#------------------------------------------Original potential factor-------------------------------------------

def potential(n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*pow(n,2))
    
    t4=  (alphaI*F/n**2) * np.exp(-3/n)
    
    
    return -t1 - t2 + t3 + t4  +I(F)/4.0


def potential_schro(n):
    
    ro=1.0/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=   (alphaI*F/n**2) * np.exp(-3/n)
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    return -t1 - t2 + t3 + t4 - t5 +I(F)/4.0


#------------------------------------------Wave number withoud dissipation--------------------------------


def kappa_C(n):

    return np.sqrt(2*mu*abs(potential_schro(n)))/hbar

def kappa(n):
    
    return np.sqrt(2*mu*abs(potential(n)))/hbar

def k():
    
    return np.sqrt(2*mu*I(F)/4)/hbar


#------------------------------------------Energy loss--------------------------------

def DE_func(T1,x,T2):
    
    RES=np.zeros(len(x))
    
    for i in range(len(x)):
        func= lambda n: gamma*kappa(n)*hbar/mu
       
        if(x[i]<=T1):
           RES[i]=0
               
        if(x[i]>T1):
            RES[i]=quad(func,T1,x[i])[0]
        
        if(x[i]>T2):
            RES[i]=quad(func,T1,T2)[0]
        
    return RES


def DE_func_C(T1_C,x,T2_C):
    
    RES_C=np.zeros(len(x))
    

    for i in range(len(x)):
        func= lambda n: gamma*kappa_C(n)*hbar/mu
        
        if(x[i]<=T1_C):
            RES_C[i]=0
        
        if(x[i]>T1_C):
            RES_C[i]=quad(func,T1_C,x[i])[0]
        
        if( x[i]>T2_C):
            RES_C[i]=quad(func,T1_C,T2_C)[0]


    return RES_C


def DE_func_num(T1,x):
    
    if(x>T1):
        func= lambda n: gamma*kappa(n)*hbar/mu
        return quad(func,T1,x)[0]
    if(x<T1):
        return 0

def DE_func_num_C(T1_C,x):
    
    if(x>T1_C):
        func= lambda n: gamma*kappa_C(n)*hbar/mu
        return quad(func,T1_C,x)[0]
    if(x<=T1_C):
        return 0

#Interpolando la funcion de dissipacion---------------------------------

x=np.linspace(0,100)
F=float(f[0])
T1=float(TP1[0])
T2=float(TP2[0])
T1_C=float(TP1_C[0])
T2_C=float(TP2_C[0])

DE=DE_func(T1,x,T2)
DE_C=DE_func_C(T1_C,x,T2_C)

interpol=interp1d(x,DE,"cubic")
interpol_C=interp1d(x,DE_C,"cubic")

plt.plot(x,DE,label="original function")
plt.plot(x,interpol(x),"k--",label="interpolation")

plt.plot(x,DE_C,label="original corrected function")
plt.plot(x,interpol_C(x),"r--",label="corrected interpolation")
plt.show()
plt.close() 

TP1[0]


#-------------------Dissipation potential--------------------------------


def disip_potential(n,inter):
    
    return potential(n)+ inter(n)

def disip_potential_C(n,inter_C):
    
    return potential_schro(n) + inter_C(n)

#---------------- Dissipative potential--------------------------

def disip_kappa(n):
    
    return np.sqrt(2*mu*(disip_potential(n)))/hbar

def disip_kappa_C(n):
    
    return np.sqrt(2*mu*(disip_potential_C(n)))/hbar

