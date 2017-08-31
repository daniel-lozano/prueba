import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from sympy.abc import r,q,a,b
from scipy.special import erf


#print(integrate((1-(a*q)**12)*exp(-(b*q)**2)*(1/q)*sin(q*r),(q,0,oo)))

a=0.316
b=0.675


c0=-np.sqrt(np.pi)*pow(a,12)/(4*pow(b,13))
c1=10395.0/32
c2=-17325/(64*pow(b,2))
c3=3465/(64*pow(b,4))
c4=-495/(128*pow(b,6))
c5=55/(512*pow(b,8))
c6=-1/(1024*pow(b,10))


r=np.linspace(0.1,0.55,1000)


alphaI=7.2
alphaN=11.08
F=0.8
m=0
Io=0.58


def potential(F,r):
    
    t1= -(alphaI*F/r**2)+F*r


    return t1


def I(F):

    return Io + (F**2)/(2*(alphaN-alphaI))






def func(r):
    return (c0*(c1+c2*r**2+c3*r**4+c4*r**6+c5*r**8+c6*r**10)*np.exp(-(r**2)/(4*b**2)) +0.5*np.pi*erf(r/(2*b))/r)

Vfs=-2*func(r)+2*potential(F,r)
V=-2/r+2*potential(F,r)

plt.plot(r,Vfs,label="finite size")
plt.plot(r,V,label="Coulomb")
plt.xlim(0,max(r))
#plt.ylim(-4,4)

plt.legend()
plt.show()


