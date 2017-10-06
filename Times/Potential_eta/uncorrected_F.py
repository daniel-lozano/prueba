import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



alphaI=7.2
alphaN=11.08
F=0.6
m=0
Io=0.58
B=0.51/2

def potential(F,n,A,B):
    
    t1= (1-(1+m)*np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4= A * np.exp(-3*B/n) * (alphaI*F/n**2) +I(F)/4.0
    
    
    return -t1 - t2 + t3 + t4




def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)


F=[0,0.2,0.4,0.6,0.8]

n=np.linspace(0.01,10,1000)
V0=potential(F[0],n,1,1)
V1=potential(F[1],n,1,1)
V2=potential(F[2],n,1,1)
V3=potential(F[3],n,1,1)
V4=potential(F[4],n,1,1)


plt.plot(n,V0,label="$ F=$"+str(F[0]) )
plt.plot(n,V1,label="$ F=$"+str(F[1]) )
plt.plot(n,V2,label="$ F=$"+str(F[2]))
plt.plot(n,V3,label="$ F=$"+str(F[3]))
plt.plot(n,V4,label="$ F=$"+str(F[4]))
plt.ylim(-1.5,1.5)

plt.title("$ Potentials $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)
plt.ylabel("$ V(\eta,F)+I_P(F)/4 $",size=15)


plt.legend(loc=1)
plt.savefig("pot_parab_F_uncorrected.png")
plt.plot()
plt.show()
