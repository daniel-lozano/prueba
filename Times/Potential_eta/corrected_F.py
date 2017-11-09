import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



alphaI=7.2
alphaN=11.08
F=0.6
m=0
Io=0.58
B=0.51/2




def potential(F,n):
    
    t1= cloud(n)#(2-(1+m)*np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2) +I(F)/4.0
    
    
    return -t1 - t2 + t3 + t4


def cloud(n):
    
    Resp=np.zeros(len(n))
    
    for i in range(len(n)):
        
        if(n[i]<= 2*B ):
            Resp[i]=(2-(1+m)*np.sqrt(I(F))/2.0)/(2*n[i])-(1/(4*B))*(1.5-(n[i]**2)/(8*B**2))
        if(n[i]> 2*B ):
            Resp[i]=(1-(1+m)*np.sqrt(I(F))/2.0)/(2*n[i])
    return Resp



def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)

def potential_schro(F,n):
    
    Z=2
    
    ro=1/(2*Z)
    
    t1= (1-(1+m)*np.sqrt(I(F))/2.0)/(2*n)#(2-(1+m)*np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2) +I(F)/4.0
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5


#F=[0,0.2,0.4,0.6,0.8]

n=np.linspace(0.01,10,1000)
V0=potential(0,n)
V1=potential(0.2,n)
V2=potential(0.4,n)
V3=potential(0.6,n)
V4=potential(0.8,n)


plt.plot(n,V0,label="$ F=$"+str(0))#str(F[0]) )
plt.plot(n,V1,label="$ F=$"+str(0.2))#str(F[1]) )
plt.plot(n,V2,label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3,label="$ F=$"+str(0.6))#str(F[3]))
plt.plot(n,V4,label="$ F=$"+str(0.8))#str(F[4]))
plt.ylim(-1.5,1.5)

plt.title("$ Potentials $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)
plt.ylabel("$ V(\eta,F) +I_P(F)/4 $",size=15)


plt.legend(loc=1)
plt.savefig("pot_parab_F_corrected.png")
plt.plot()
plt.show()



n=np.linspace(0.01,10,1000)
V0=potential_schro(0,n)
V1=potential_schro(0.2,n)
V2=potential_schro(0.4,n)
V3=potential_schro(0.6,n)
V4=potential_schro(0.8,n)


plt.plot(n,V0,label="$ F=$"+str(0))#str(F[0]) )
plt.plot(n,V1,label="$ F=$"+str(0.2))#str(F[1]) )
plt.plot(n,V2,label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3,label="$ F=$"+str(0.6))#str(F[3]))
plt.plot(n,V4,label="$ F=$"+str(0.8))#str(F[4]))
plt.ylim(-1.5,1.5)

plt.title("$ Potentials $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)
plt.ylabel("$ V(\eta,F)+I_P(F)/4 $",size=15)


plt.legend(loc=1)
plt.savefig("pot_parab_F_corrected_schro.png")
plt.plot()
plt.show()
