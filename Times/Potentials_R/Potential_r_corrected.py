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

alphaI=7.2
alphaN=11.08
F=0.08
m=0
Io=0.58



"""
-----------------------------FUNCTIONS-------------------------------
"""
def I(F):
    
    return Io + (F**2)/(2*(alphaN-alphaI))

def nucleus_func(r,a,b):
    
    c0=-pow(a,12)/(8*pow(b,13)*pow(np.pi,1.5))
    c1=10395.0/32
    c2=-17325/(64*pow(b,2))
    c3=3465/(64*pow(b,4))
    c4=-495/(128*pow(b,6))
    c5=55/(512*pow(b,8))
    c6=-1.0/(1024*pow(b,10))
    
    Norm=1/0.0795#(np.pi)**2
    
    return (c0*(c1+c2*r**2+c3*r**4+c4*r**6+c5*r**8+c6*r**10)*np.exp(-(r**2)/(4*b**2)) +erf(r/(2*b))/(4*np.pi*r))*Norm

def e_func(r,leng):
    
    
    V=np.zeros(leng)
    
    for i in range(leng):
        
        if(r[i]<=A):
            V[i]=0#1/B +( (B**2-A**2)/2 -(A**3)*(1/A-1/B))
        if(r[i]<=B and r[i]>A):
            V[i]=+1/B +( (B**2-r[i]**2)/2 -(A**3)*(1/r[i]-1/B))
        if(r[i]>B):
            V[i]=1/r[i]
    
    return V

def e_func_sch(r):
    Z=2
    ro=1/(2*Z)
    return -1/r +(1/r)*(1+r/(2*ro))*np.exp(-r/ro)


"""
-----------------------------Graphing potentials-------------------------------
"""
a=0.316E-5
b=0.675E-5
r=np.linspace(0.001,8,1000)
B=(0.529/2)#*1E5
A=0#2E-5
F=-1

factor=1#1.0E-5 #added to account for the units used in the AU

#Vfs=-1/r#-2*nucleus_func(r,a,b) + e_func(r,len(r))+ F*r*factor -(alphaI*F/(factor*r)**2)

V1=-2/r+e_func(r,len(r))
V2= -2/r+e_func(r,len(r)) +F*r
V3=-2/r+e_func(r,len(r)) + F*r -(alphaI*F/(r)**2)#*np.exp(-3/r)
V4=-2/r+e_func(r,len(r)) + F*r -(alphaI*F/(r)**2)*np.exp(-3/r)



plt.plot(r,V1,label="$ -2/r +V_{cloud} $")
plt.plot(r,V2,label="$ -2/r +V_{cloud} + r F_{p} $")
plt.plot(r,V3,label="$ -2/r +V_{cloud} + r F_{p} + V_{pol}$")
plt.plot(r,V4,"k",label="$ -2/r +V_{cloud} + r F_{p} + V_{pol}exp(-3/r) $")


plt.title("$ Effect\ of\ corrections\ in\ potential $")
plt.xlim(0,max(r))
plt.ylim(-20,10)
plt.xlabel("$ r (Angstroms)  $",size=15)
plt.ylabel("$ V(r) $",size=15)


plt.legend(loc=4)

plt.savefig("Effect_of_polirazation_corrected.png")
plt.show()
plt.close()

"""
-----------------------------Graphing potentials schro-------------------------------
"""
a=0.316E-5
b=0.675E-5
r=np.linspace(0.001,8,1000)
B=(0.529/2)#*1E5
A=0#2E-5
F=-1

factor=1#1.0E-5 #added to account for the units used in the AU

#Vfs=-1/r#-2*nucleus_func(r,a,b) + e_func(r,len(r))+ F*r*factor -(alphaI*F/(factor*r)**2)


V3=-1/r+ F*r -(alphaI*F/(r)**2)*np.exp(-3/r)
V4=-2/r+e_func_sch(r) + F*r -(alphaI*F/(r)**2)*np.exp(-3/r)


plt.plot(r,V3,label="$ -1/r + r F_{p} + V_{pol}exp(-3/r)$")
plt.plot(r,V4,"k",label="$ -2/r +V_{schr} + r F_{p} + V_{pol}exp(-3/r) $")


plt.title("$ Effect\ of\ corrections\ in\ potential\ Schro $")
plt.xlim(0,max(r))
plt.ylim(-20,10)
plt.xlabel("$ r (Angstroms)  $",size=15)
plt.ylabel("$ V(r) $",size=15)


plt.legend(loc=4)

plt.savefig("Effect_of_polirazation_corrected_schro.png")
plt.show()
plt.close()









