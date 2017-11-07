import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy import integrate


alphaI=7.2
alphaN=11.08
m=0
Io=0.58
B=0.51/2

F=0

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)



def potential_schro(n):
    
    Z=2
    
    ro=1/(2*Z)
    
    t1= (1-(1+m)*np.sqrt(I(F))/2.0)/(2*n)#(2-(1+m)*np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2) +I(F)/4.0
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5
F=0.01
n=np.linspace(0.1,60,1000)
plt.plot(n,potential_schro(n))
plt.ylim(-0.5,0.5)
plt.show()




x=np.linspace(0.04,0.12,20)

for i in x:
    F=i

    print("Value of F=", F)
    print("[",brentq(potential_schro,0.1,4),",",brentq(potential_schro,4,30),"]")











