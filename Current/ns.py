import matplotlib.pyplot as plt
import numpy as np


T=np.linspace(0,9.25) # K
mc2=2*0.51E6 #eV
lam=47.0E-3 #nm
TD=276.0 #k
N0=19.87 # 1/eV
g=1/(N0*np.log(1.13*TD/9.25))    #eV
e=2.0#*1.6E-19
print("g=",g)


ns=(mc2/(4*e*np.pi*pow(lam,2)))*(1-pow(T,4)*np.exp(4.0/(g*N0))/pow(1.13*TD,4))

plt.plot(T,ns)
plt.xlabel("$ T $")
plt.ylabel("$ \eta_s(T)  $")
plt.xticks([9.25],["$ T_c $"])
plt.ylim(0,max(ns)*1.05)
plt.savefig("ns.png")
plt.show()
