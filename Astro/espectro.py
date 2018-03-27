import numpy as np
import matplotlib.pyplot as plt
from sys import argv

Labels=["N_2","H_2","O_2","CO","NO"]

R=[109.4,74.16,120.7,112.8,115.1]
Be=[2.01,60.8,1.45,1.931,1.705]
De=[5.8E-6,1.6E-2,4.8E-6,6E-6,0.5E-6]
alpha_e=[0.017,3.06,0.016,0.017,0.017]
beta_e=alpha_e#[0.0,0.0,0.0,0.0,0.0]
w_e=[2359.0,4401,1580.0,2170.0,1904.0]
wxe=[14.3,121.3,12.0,13.29,14.08]

mole=int(argv[1])

#-------------------------------------------Definiendo funciones--------------------------------------------

def Energy(n,l,mol):
    
    return w_e[mol]*(n+0.5) - wxe[mol]*(n+0.5)**2 + Be[mol]*l*(l+1) + De[mol]*pow(l*(l+1),2) -(n+0.5)*( alpha_e[mol]*l*(l+1)+beta_e[mol]* pow(l*(l+1),2) )

def delta(n1,n2,l1,l2,mol):

    return Energy(n2,l2,mol)-Energy(n1,l1,mol)

print delta(1,2,1,2,mole)



line=[0,1.0]
ones=np.ones(2)

fig=plt.figure(figsize=(15,5))
ax=fig.add_subplot(111)

n_plots=0

for n in range(3): #movimiento de n
    print n
    for dn in range(3): #cambio de n
    
        for l in range(3):  #   moviemiento de l
            
            p0=delta(n,n+dn,l,l+0,mole)#
            p1=delta(n,n+dn,l,l+1,mole)#
            p2=delta(n,n+dn,l,l+2,mole)#
            
            if(p0!=0):
                n_plots+=1
                P0=ones*p0
                ax.plot(P0,line)
            
            n_plots+=1
            P1=ones*p0
            ax.plot(P1,line)

            n_plots+=1
            P2=ones*p0
            ax.plot(P2,line)


print "number of plots=",n_plots

plt.xlim(2000,7000)
plt.title("$ "+Labels[mole]+" $",size=20)
plt.xlabel("$  \mu m $", size=15)
plt.yticks([])
plt.show()



