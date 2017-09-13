import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

kb=8.6E-5
Tc=9.25
Do=30.5E-4 # eV

Tl=np.linspace(0.1*Tc,0.6*Tc)
Th=np.linspace(0.6*Tc,Tc)
T=np.linspace(0.1*Tc,Tc)


def high(T):
    return  3.07*kb*Tc*np.sqrt(1-T/Tc)

def low(T):

    return Do*(1-(2*np.pi*kb*T*np.exp(-Do/(kb*T)))/Do)

plt.plot(T,low(T))
plt.plot(T,high(T))
plt.show()
