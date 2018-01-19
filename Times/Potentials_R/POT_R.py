import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from sympy.abc import r,q,a,b
from scipy.special import erf
import sys


#print(integrate((1-(a*q)**12)*exp(-(b*q)**2)*(1/q)*sin(q*r),(q,0,oo)))
"""
---------------------------CONSTANTS AND VARIABLES------------------------------
"""
B=0.529/2
A=2E-5

a=0.316#0
b=0.675#0

r=np.linspace(1E-20,15,10000)

alphaI=0.28125
alphaN=1.38
m=0
Io=0.904
B=0.51/2
Z=2
F=0.08


"""
-----------------------------FUNCTIONS-------------------------------
"""
def I(F):
    
    return Io + (F**2)/(2*(alphaN-alphaI))

def e_func(r,F,theta):
    Z=2.0
    ro=1/(2*Z)
    return -(Z-1)/r  + F*r*np.cos(theta) -(np.cos(theta)*alphaI*F/(r)**2)*np.exp(-3/r) +I(F)

def e_func_sch(r,F,theta):
    Z=2.0
    ro=1/(2*Z)
    return -(Z-1)/r -(1/r)*(1+r/(2*ro))*np.exp(-r/ro)  + F*r*np.cos(theta) -(np.cos(theta)*alphaI*F/(r)**2)*np.exp(-3/r) +I(F)



"""
-----------------------------Graphing potentials schro-------------------------------
"""

r=np.linspace(0.01,100,100)
F1=-0.04
F2=-0.11


Theta=[0,30,60]
for i in range(len(Theta)):
    theta=Theta[i]
    E4=np.ones(len(r))*I(F1)*-1
    plt.plot(r,E4,"g",label="F="+str(F1))
    E11=np.ones(len(r))*I(F2)*-1
    plt.plot(r,E11,"b",label="F="+str(F2))
    
    
    V4=e_func(r,F1,theta*np.pi/180.0)
    V11=e_func(r,F2,theta*np.pi/180.0)
    
    VC4=e_func_sch(r,F1,theta*np.pi/180.0)
    VC11=e_func_sch(r,F2,theta*np.pi/180.0)
    print("Ultima posicion",VC4[-1])
    
    plt.plot(r,V4,"k",label="0.04 Uncorrected")
    plt.plot(r,V11,"r",label="0.11 Uncorrected")
    plt.plot(r,VC4,"k--",label="0.04 Corrected")
    plt.plot(r,VC11,"r--",label="0.11 Corrected")
    plt.title(" theta= "+str(theta))
    
    
    
    
    
    plt.legend()
    plt.plot()

    plt.ylim(-5,5)
    plt.xlabel("$ r (Angstroms)  $",size=15)
    plt.ylabel("$ V(r) $",size=15)
    plt.legend(loc=4)
    #plt.savefig("r_plot_theta="+str(theta)+".png")
    plt.show()
    plt.close()



