import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



phi=np.linspace(0,50*np.pi,1000)

beta=1.1
E=0.5
ro=1

BETA=[0.1,0.25,0.5,0.9,0.98,0.99,0.999,1,1.001,1.01,1.02,1.03,1.04,1.05,np.pi/3,np.pi/2]

BETA.sort()


print(BETA)
for i in BETA:


    beta=i
    theta=beta*phi
 

    u=1+E*np.sin(beta*phi)#se puede tener soluciones que tengan + o - dependiendo
    r=ro/(u)

    ax=plt.subplot(111,projection='polar')
    ax.plot(phi,r)
    ax.set_rmax(3)#(max(r)-1)
    ax.set_title("$ Orbita, $ B=  "+str(beta),va="bottom")
    
    plt.show()


