import numpy as np
import matplotlib.pyplot as plt


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
