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


print(c0*c1,c0*c2,c0*c3,c0*c4,c0*c5,c0*c6)

r=np.linspace(0.1,15,1000)


def func(r):
    return (c0*(c1+c2*r**2+c3*r**4+c4*r**6+c5*r**8+c6*r**10)*np.exp(-(r**2)/(4*b**2)) +0.5*np.pi*erf(r/(2*b))/r)

V=-2*func(r)

plt.plot(r,V,label="finite size")
plt.plot(r,-1/r,label="Coulomb")
plt.xlim(0,max(r))
plt.ylim(-4,4)

plt.legend()
plt.show()


