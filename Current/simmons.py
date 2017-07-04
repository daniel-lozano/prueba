import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate




a=20.7989  # Angstrom
Vo=1.8 +3.4E-4# eV Vo+Delta
m=0.51E6     # eV
hbar=1973.27 # eV*Angstrom =[hbar*c]
h=2*np.pi*hbar # eV*Angstrom
hs=4.13E-15
e=1.6E-19 # Coulomb
kb=8.6E-5# eV/Kelvin
T=1.5#4 Kelvin
beta=1.0/(kb*T)
kt=1.0/beta

print("A=",a, "Angstroms")
print("Vo+Deltao=",Vo, "eV")
print("m=",m,"eV")
print("hbar=",hbar,"eV*Angstrom")
print("e=",e,"Coulomb")
print("beta=",beta, "1/eV")





mu=0#chemical potential, no particles are added

g=kb*T*np.log(3) #groundstate taken at 1/2 of Bose-Einsteins distribution

EV=np.linspace(-1,1,1000)#(-50E15,10E15,1000)*2*e#applied potential in eV

V=EV/(2*e) #Applied potential in V


A=(4*np.pi*a*np.sqrt(2*m))/h
B=A/(2*np.sqrt(Vo))



J=e*1E16*( Vo*np.exp(-A*np.sqrt(Vo))- (Vo+EV)* np.exp(-A*np.sqrt(Vo+EV)) )/(2*np.pi*hs*pow(a,2))



plt.plot(EV,J)
plt.show()
plt.close()
