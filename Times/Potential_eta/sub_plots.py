import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



alphaI=0.28125
alphaN=1.38
m=0
Io=0.904
B=0.51/2
Z=2.0
F=0


def potential(F,n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2.0)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4= np.exp(-3/n) * (alphaI*F/n**2) #+I(F)/4.0
    
    
    return -t1 - t2 + t3 + t4


def potential_C(F,n):
    
    t1= cloud(n)#(2-(1+m)*np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2) #+I(F)/4.0
    
    
    return -t1 - t2 + t3 + t4


def cloud(n):
    
    Resp=np.zeros(len(n))
    
    for i in range(len(n)):
        
        if(n[i]<= 2*B ):
            Resp[i]=(2-(1+m)*np.sqrt(I(F))/2.0)/(2*n[i])-(1/(4*B))*(1.5-(n[i]**2)/(8*B**2))
        if(n[i]> 2*B ):
            Resp[i]=(1-(1+m)*np.sqrt(I(F)/2.0))/(2*n[i])
    return Resp



def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/2

def potential_schro(F,n):
    
    
    
    ro=1.0/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2.0)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2) #+I(F)/4.0
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5




#F=[0,0.2,0.4,0.6,0.8]

plt.figure(figsize=(10,5))
n=np.linspace(0.01,5,1000)

V0=potential(0,n)
V1=potential(0.2,n)
V2=potential(0.4,n)
V3=potential(0.6,n)
V4=potential(0.8,n)

V0C=potential_C(0,n)
V1C=potential_C(0.2,n)
V2C=potential_C(0.4,n)
V3C=potential_C(0.6,n)
V4C=potential_C(0.8,n)

plt.subplot(1,1,1)
plt.grid(True)


plt.plot(n,V1,"b",label="$ F=$"+str(0.2))#str(F[1]) )
plt.plot(n,V2,"r",label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3,"g",label="$ F=$"+str(0.6))#str(F[3]))
plt.plot(n,V4,"c",label="$ F=$"+str(0.8))#str(F[3]))




plt.plot(n,V1C,"b--")#str(F[1]) )
plt.plot(n,V2C,"r--")#str(F[2]))
plt.plot(n,V3C,"g--")#str(F[3]))
plt.plot(n,V4C,"c--")#str(F[3]))

plt.ylabel("$ V(\eta,F) $",size=15)
plt.ylim(-1,0.5)
plt.legend(loc=2)
plt.title("$ Comparing\ potentials $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)



plt.legend(loc=1)
plt.savefig("subplot.png")
plt.plot()
plt.show()

"""
------------------ POTENTIAL WITH SCHRODINGER----------------------------------------------
"""
plt.figure(figsize=(15,10))
n=np.linspace(0.01,100,1000)
V1=potential(0.04,n)
V1C=potential_schro(0.04,n)
plt.plot(n,V1C,"r",label="$ F=$"+str(0.04))#str(F[1]))
plt.plot(n,V1,"k",label="$ F=$"+str(0.04)+"Uncorrected")#str(F[1]))
plt.title("F=0.04")
plt.ylim(-0.4,0)
plt.xlim(0,5)
plt.legend()
plt.show()

"""
------------------ POTENTIAL WITH SCHRODINGER----------------------------------------------
    """


plt.figure(figsize=(15,10))
n=np.linspace(0.01,10,100)

V2=potential(0.4,n)
V3=potential(0.6,n)
V4=potential(0.8,n)

plt.subplot(1,1,1)
plt.grid(True)


ONES=np.ones(len(n))*(-I(0.4)/4)
plt.plot(n,ONES,"k.")
plt.plot(n,V2,"k",label="$ F=$"+str(0.4))#str(F[2]))

ONES=np.ones(len(n))*(-I(0.6)/4)
plt.plot(n,ONES,"b.")
plt.plot(n,V3,"b",label="$ F=$"+str(0.6))#str(F[3]))

ONES=np.ones(len(n))*(-I(0.8)/4)#
plt.plot(n,ONES,"g.")
plt.plot(n,V4,"g",label="$ F=$"+str(0.8))#str(F[3]))

V2C=potential_schro(0.4,n)
V3C=potential_schro(0.6,n)
V4C=potential_schro(0.8,n)

plt.plot(n,V2C,"k--")#,label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3C,"b--")#,label="$ F=$"+str(0.6))#str(F[3]))
plt.plot(n,V4C,"g--")#,label="$ F=$"+str(0.8))#str(F[3]))



plt.ylim(-0.5,-0.15)
plt.title("$ Hydrogen\ atom\ solution\ $",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)
plt.ylabel("$ V(\eta,F) $",size=15)


plt.legend(loc=1)
plt.savefig("subplot_schro.png")
plt.plot()
plt.show()

"""
    ------------------ POTENTIALS COMPARISON----------------------------------------------
    """

plt.figure(figsize=(10,5))
n=np.linspace(0.01,10,100)

V2=potential(0.4,n)
V3=potential(0.8,n)

V2C=potential_C(0.4,n)
V3C=potential_C(0.8,n)

V2S=potential_schro(0.4,n)
V3S=potential_schro(0.8,n)


plt.subplot(1,1,1)
plt.grid(True)

plt.plot(n,V2,"r",label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3,"k",label="$ F=$"+str(0.8))#str(F[3]))

plt.plot(n,V2C,"r--")#,label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3C,"k--")#,label="$ F=$"+str(0.6))#str(F[3]))


plt.plot(n,V2S,"r-.",linewidth=3)#,label="$ F=$"+str(0.4))#str(F[2]))
plt.plot(n,V3S,"k-.",linewidth=3)#,label="$ F=$"+str(0.6))#str(F[3]))


plt.ylim(-0.5,-0.15)
plt.title(" Cloud Potential vs. Hydrogen atom Potential ",size=20)
plt.xlabel("$ \eta\ (a.u.) $",size=15)
plt.ylabel("$ V(\eta,F) $",size=15)
plt.legend(loc=4)

plt.savefig("subplot_comparing.png")
plt.plot()
plt.show()






