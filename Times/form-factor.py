import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from sympy import *
from sympy.abc import r,q
from sympy.integrals.transforms import fourier_transform



#a=0.316 # fermi
#b=0.675 # fermi





Integ=integrate((1-(a*q)**12)*exp(-(b*q)**2)*q*sin(q*r),(q,0,oo))

value=integrate(q*sin(r*q),(q,0,oo))

print("value=",Integ)


#a=0.316 # fermi
#b=0.675 # fermi

def func(r):

    return (1/(np.sqrt(np.pi)**3))*  ((-4.5E-5)*(0.027*r**12 - 1.939*r**10 + 48.606*r**8 - 531.516*r**6 + 2542.803*r**4 - 4634.259*r**2  + 135135.0/64)*np.exp(-0.5487*r**2) + 0.4065* np.exp(-0.5487*r**2))


r=np.linspace(0.1,15)
V=-func(r)



gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])
ax3 = plt.subplot(gs[1,: ])

ax1.set_title("finite size")
ax1.plot(r,V,"r")
ax1.set_ylabel("$ V(ev)   $")
ax1.set_xlabel("$ r(rm)   $")


ax2.plot(r,-1/r)
ax2.set_title("Coulomb")
ax2.set_xlabel("$ r(rm)   $")



plt.plot(r,-1/r)
ax3.plot(r,V,"r")
ax3.set_ylabel("$ V(ev)   $")
ax3.set_xlabel("$ r(rm)   $")


plt.savefig("finite_vs_coulomb.png")
plt.show()
plt.close()
