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
    
    t1= (1-(1+m)*np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4= np.exp(-3/n) * (alphaI*F/n**2) +I(F)/4.0
    
    
    return -t1 - t2 + t3 + t4



def potential_C(F,n):
    
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


#F=[0,0.2,0.4,0.6,0.8]

plt.figure(figsize=(10,5))
n=np.linspace(0.01,5,1000)
V0=potential(0,n)
V1=potential(0.2,n)
V2=potential(0.4,n)
V3=potential(0.6,n)
V4=potential(0.8,n)

plt.subplot(1,2,1)
plt.grid(True)

#plt.plot(n,V0,"k",label="$ F=$"+str(0))#str(F[0]) )
plt.plot(n,V1,"b",label="$ F=$"+str(0.2))#str(F[1]) )
plt.plot(n,V2,"r",label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3,"g",label="$ F=$"+str(0.6))#str(F[3]))
plt.plot(n,V4,"c",label="$ F=$"+str(0.6))#str(F[3]))
plt.ylim(-0.3,0.5)


plt.title("$ Uncorrected\ potential $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)
plt.ylabel("$ V(\eta,F)+I_P(F)/4 $",size=15)
plt.legend(loc=4)







n=np.linspace(0.01,5,1000)
V0=potential_C(0,n)
V1=potential_C(0.2,n)
V2=potential_C(0.4,n)
V3=potential_C(0.6,n)
V4=potential_C(0.8,n)

plt.subplot(1,2,2)
plt.grid(True)

#plt.plot(n,V0,"k--",label="$ F=$"+str(0))#str(F[0]) )
plt.plot(n,V1,"b",label="$ F=$"+str(0.2))#str(F[1]) )
plt.plot(n,V2,"r",label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3,"g",label="$ F=$"+str(0.6))#str(F[3]))
plt.plot(n,V4,"c",label="$ F=$"+str(0.6))#str(F[3]))


plt.ylim(-0.3,0.5)

plt.title("$ Corrected\ potential $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)



plt.legend(loc=4)
plt.savefig("subplot.png")
plt.plot()
plt.show()
