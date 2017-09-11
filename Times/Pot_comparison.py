import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from sympy.abc import r,q,a,b
from scipy.special import erf
import sys


#print(integrate((1-(a*q)**12)*exp(-(b*q)**2)*(1/q)*sin(q*r),(q,0,oo)))

a=0.316E-5
b=0.675E-5

B=0.529/2
A=1.75E-5


c0=-np.sqrt(np.pi)*pow(a,12)/(4*pow(b,13))
c1=10395.0/32
c2=-17325/(64*pow(b,2))
c3=3465/(64*pow(b,4))
c4=-495/(128*pow(b,6))
c5=55/(512*pow(b,8))
c6=-1/(1024*pow(b,10))


r=np.linspace(1E-6,2,10000)

alphaI=7.2
alphaN=11.08
F=0.08
m=0
Io=0.58



def I(F):

    return Io + (F**2)/(2*(alphaN-alphaI))

def nucleus_func(r,a,b):
    
    return (c0*(c1+c2*r**2+c3*r**4+c4*r**6+c5*r**8+c6*r**10)*np.exp(-(r**2)/(4*b**2)) +0.5*np.pi*erf(r/(2*b))/r)

def e_func(r,leng):


    V=np.zeros(leng)

    for i in range(leng):
        
        if(r[i]<=A):
            V[i]=0
        if(r[i]<=0.5 and r[i]>A):
            V[i]=-1/B -( (B**2-r[i]**2)/2 -(A**3)*(1/r[i]-1/B))
        if(r[i]>B):
            V[i]=-1/r[i]

    return V







Vfs=-2*nucleus_func(r,a,b) #+ e_func(r,len(r))

V=-2/r

plt.plot(r,Vfs-V,label="finite size")
plt.plot(r,V,label="Coulomb")
plt.xlim(min(r),max(r))
plt.ylim(-15,4)
plt.xlabel("$ r (Angstroms)  $",size=15)
plt.ylabel("$ V(r) $",size=15)
plt.title("$ V_{FS}\ Vs.\ Coulomb  $")
plt.savefig("comparing_potentials.png")
plt.legend(loc=4)
plt.show()
plt.close()

Vfs=-2*nucleus_func(r,a,b) + e_func(r,len(r))

V=-2/r

plt.plot(r,Vfs-V,label="finite size")
plt.plot(r,V,label="Coulomb")
plt.xlim(min(r),max(r))
plt.ylim(-10,0)
plt.xlabel("$ r (Angstroms)  $",size=15)
plt.ylabel("$ V(r) $",size=15)
plt.title("$ V_{FS}+V_e\ Vs.\ Coulomb  $")
plt.savefig("comparing_potentials_final.png")
plt.legend(loc=4)
plt.show()
plt.close()



"""
Vfs=-2*I(F)*func(r)-(alphaI*F/r**2)+F*r
V=-2*I(F)/abs(r1) -(alphaI*F/r1**2)+F*r1

plt.plot(r,Vfs,label="finite size")
plt.plot(r1,V,label="Coulomb")
#plt.xlim(0,max(r))
plt.ylim(-4,4)
plt.savefig("comparing_potentials.png")
plt.legend()
plt.show()
"""

