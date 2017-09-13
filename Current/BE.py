import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

kb=8.6E-5


T=[0.5, 1, 2, 5]


def func(E,T):

    return 1/(np.exp(E/(kb*T))-1)

E=np.linspace(0.3*kb,15*kb,1000)

A=np.ones(len(E))*0.5
B=[0,0.5]
xtl=[]
xt=[]



for i in range(len(T)):
    

    f=func(E,T[i])
    l="$ T_"+str(i)+"= $"
    plt.plot(E,f,label=l+str(T[i])+"$ \ K $ ")
    
    C=np.ones(2)*kb*T[i]*np.log(3)


#labeling ticks

    xtl.append("$ \eta_"+str(i)+"$")
    xt.append(kb*T[i]*np.log(3))
    plt.plot(C,B,"k--")


plt.legend()
plt.plot(E,A,"k--")
plt.title("B-E distribution")
plt.xlim(min(E)*0.8,max(E)/3)
plt.ylim(-0.3,max(func(E,T[-1])))
plt.yticks([0.5],["1/2"])
plt.xticks(xt,xtl,rotation="horizontal",size=15)

plt.savefig("BE_distribution.png")

plt.show()
plt.close()

T=np.linspace(0,1,100)
D=3.07*(15.5E-4)*np.sqrt(1-T)


plt.plot(T,D)
plt.xlim(0,1.2)
plt.ylim(0,max(D)+max(D)*0.2)
plt.yticks([D[0]],[" $ \Delta_0 $ "],size=15)
plt.xticks([1],[" $ T/T_c  $"],size=15)
plt.title("$ \Delta(T) $")
plt.savefig("Delta.png")
plt.show()

g=1E-5 #eV/e


N=D/(2*g)+0.5
plt.plot(T,N)
plt.show()


















