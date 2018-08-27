import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad

from sys import argv



alphaI=0.28125
alphaN=1.38
m=0
Io=0.904#
B=0.51/2#0.51/2

mu=1.0 #a.u.
hbar=1.0 #a.u.

h=2.0*np.pi
c=137.0
lam=735.0/0.0529

Z=2
omega=h*c/lam #
number=0.0#float(argv[1])
#gamma=float(argv[1])
#print("gamma_D=",gamma)


f=np.linspace(0.04,0.11,7)
angle=np.zeros(len(f))

Factor=24.18884 #time[as]/a.u.

#------------------------------------------Energy factor-------------------------------------------
#D=0 #Dissipation

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)


def V(r,r_0,cte):
    return -1/r-(1/r)*(1+r/(2*r_0))*np.exp(-r/r_0)*cte
r=np.linspace(B,1E4)

#
#R=np.linspace(0.01,B)
#for i in range(len(f)-1):
#    func_r=2*B*(pow(R,4)-pow(B*R,2) -pow(R,4)*V(R,B)/I(f[i]) )**(-0.5)
#    plt.plot(R,func_r,label="$ F =$"+str(f[i]))
#plt.legend()
#plt.show()

for i in range(len(f)):

    func=lambda r: 2*B*(r**4-(B*r)**2 -(r**4)*V(r,B,0)/I(f[i]) )**(-0.5)#*
    angle[i]=quad(func,B,np.inf)[0]#-np.pi/2

plt.plot(f,angle,label="$ \\theta_m $")
plt.title("$ v(r)=\\frac{1}{r}  $")
plt.legend()
plt.savefig("angle.png")
plt.show()




