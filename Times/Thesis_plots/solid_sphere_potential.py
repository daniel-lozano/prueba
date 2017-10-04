import numpy as np
import matplotlib.pyplot as plt

'''
a=1
z=2
r=np.linspace(a+0.1,a*5,1000)

Vss=-a/(2*(r**2-a**2))-(z-a/(2*r))/r
Vc=-2/r




plt.plot(r,Vss,label="Solid sphere")
plt.plot(r,Vc,label="coulomb")
plt.savefig("Solid_sphere_coulomb.png")
plt.legend()
plt.show()

'''

x=np.linspace(-5,5,3000)
vx=-1/abs(x)+x*15*np.cos(3*np.pi/4)
y=np.ones(len(x))*-10

plt.plot(x,vx)
plt.plot(x,y,"r--")
plt.title("$Attoclock\ potential$",size=20)
plt.xlabel("$ x $",size=15)
plt.ylabel("$ V(x) $",size=15)
plt.xticks([],[])
plt.yticks([],[])
plt.xlim(-2,2)
plt.ylim(-20,20)
plt.savefig("x_potential.png")
plt.legend()
plt.show()
