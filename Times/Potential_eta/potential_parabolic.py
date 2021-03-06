import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



alphaI=7.2
alphaN=11.08
F=0.6
m=0
Io=0.58


def potential(F,n,A,B):
    
    t1= (1-np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4= A * np.exp(-3*B/n) * (alphaI*F/n**2) +I(F)/4.0


    return -t1 - t2 + t3 + t4


def I(F):

    return Io + ((alphaN-alphaI)*F**2)/(2)




n=np.linspace(0.1,20,1000)
V0=potential(0,n,0,0)
V1=potential(F,n,0,0)
V2=potential(F,n,1,0)
V3=potential(F,n,1,1)


plt.plot(n,V0,label="$ F=0\ V_{pol}=0\ $" )
plt.plot(n,V1,label="$ V_{pol}=0\ $" )
plt.plot(n,V2,label=" $ V_{pol}\ $ ")
plt.plot(n,V3,label="$  V_{pol} \cdot exp(-3/\sigma) $ ")
plt.ylim(-1.5,1.5)

plt.title("$ Potentials $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)
plt.ylabel("$ V(\eta,F)+I_p(F)/4 $",size=15)


plt.legend(loc=1)
plt.savefig("potential_parabolic.png")
plt.plot()
plt.show()









