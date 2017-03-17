import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



phi=np.linspace(0,50*np.pi,1000)

beta=1.1
E=1/2.0
ro=1



#-----------RACIONALES---------------------------

BETA=[0.1,0.25,0.5,0.9,1,1.1,1.25,1.5,3]#[0.1,0.25,0.5,0.9,0.98,0.99,0.999,1,1.001,1.01,1.02,1.03,1.04,1.05,np.pi/3,np.pi/2]

BETA.sort()

fig = plt.figure(figsize=(10, 10))
print(BETA)
for i in range(len(BETA)):


    beta=BETA[i]
    theta=beta*phi
 

    u=1+E*np.cos(beta*phi)#se puede tener soluciones que tengan + o - dependiendo
    r=ro/(u)
    posicion="33"+str(i+1)
    ax=plt.subplot(posicion,projection='polar')
    ax.plot(phi,r)
    ax.set_rmax(3)
    ax.set_title("$  B=  $ "+str(beta),va="bottom")
    ax.set_rticks([0, 1, 2])
    ax.set_xticks([0])
plt.savefig("racionales.png")
plt.show()

plt.close()
#-----------CERCANOS A 1---------------------------

BETA=[0.9,0.95,0.98,0.99,1,1.01,1.02,1.03,1.04]#[0.1,0.25,0.5,0.9,0.98,0.99,0.999,1,1.001,1.01,1.02,1.03,1.04,1.05,np.pi/3,np.pi/2]

BETA.sort()

fig = plt.figure(figsize=(10, 10))
print(BETA)
for i in range(len(BETA)):
    
    
    beta=BETA[i]
    theta=beta*phi
    
    
    u=1+E*np.cos(beta*phi)#se puede tener soluciones que tengan + o - dependiendo
    r=ro/(u)
    posicion="33"+str(i+1)
    ax=plt.subplot(posicion,projection='polar')
    ax.plot(phi,r)
    ax.set_rmax(3)
    ax.set_title("$  B= $ "+str(beta),va="bottom")
    ax.set_rticks([0, 1, 2])
    ax.set_xticks([0])

plt.savefig("cerca_1.png")
plt.show()

plt.close()


#-----------IRRACIONALES---------------------------

BETA=[np.e/2,np.pi/3,np.pi/4,np.e/3]
labels=["$ B=  \pi/4 $","$ B= e/3 $","$ B= \pi/3 $","$ B= e/2 $"]
BETA.sort()

fig = plt.figure(figsize=(10, 10))
print(BETA)
for i in range(len(BETA)):
    
    
    beta=BETA[i]
    theta=beta*phi
    
    
    u=1+E*np.cos(beta*phi)#se puede tener soluciones que tengan + o - dependiendo
    r=ro/(u)
    posicion="22"+str(i+1)
    ax=plt.subplot(posicion,projection='polar')
    ax.plot(phi,r)
    ax.set_rmax(3)
    ax.set_title(labels[i],va="bottom")
    ax.set_rticks([0, 1, 2])
    ax.set_xticks([0])

plt.savefig("irracionales.png")

plt.show()

plt.close()













