import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad

from sys import argv



alphaI=0.28125
alphaN=1.38
m=0
Io=0.904#
B=0.51/2.0#0.51/2

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
anglec=np.zeros(len(f))


Factor=24.18884 #time[as]/a.u.

#------------------------------------------Energy factor-------------------------------------------
#D=0 #Dissipation

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)

def V(r):
    return -1./r


def Vc(r,r_0):
    return -1./r-(1./r)*(1+r/(2*r_0))*np.exp(-r/r_0)


r=np.linspace(B,1E4)

#
#R=np.linspace(0.01,B)
#for i in range(len(f)-1):
#    func_r=2*B*(pow(R,4)-pow(B*R,2) -pow(R,4)*V(R,B)/I(f[i]) )**(-0.5)
#    plt.plot(R,func_r,label="$ F =$"+str(f[i]))
#plt.legend()
#plt.show()

for i in range(len(f)):

    func1=lambda r: 2*B*(r**4-(B*r)**2 -(r**4)*V(r)/I(f[i]))**(-0.5)#*
    angle[i]=quad(func1,B,np.inf)[0]#-np.pi/2

    func2=lambda r: 2*B*(r**4-(B*r)**2 -(r**4)*Vc(r,B)/I(f[i]) )**(-0.5)#*
    anglec[i]=quad(func2,B,np.inf)[0]#-np.pi/2

plt.figure(figsize=(16,5))
plt.title("$ \\theta_{offset} $")
plt.subplot(131)
plt.plot(f,angle,label="$ -\\frac{1}{r}$")
plt.legend()
plt.subplot(132)
plt.plot(f,anglec,label="$ -\\frac{1}{r}\ +\ \mathrm{terms} $")
plt.legend()
plt.subplot(133)
plt.plot(f,(angle-anglec),label="$ \mathrm{dif} $")#*180/np.pi
plt.legend()
plt.savefig("angle.png")
plt.show()




