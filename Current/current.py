import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate




dx=20.7989  # Angstrom
e=1.6E-19 # Coulomb
m=0.51E6   # eV


Vo=(1.8 +3.4E-4)# eV, Vo+Delta

hbar=1973.27# eV*Angstrom =[hbar*c]

h=2*np.pi*hbar # eV*Angstrom
hs=6.5E-16 # eV s

kb=8.6E-5# eV/Kelvin/e
T=1.5#4 Kelvin
beta=1.0/(kb*T)

print("A=",dx, "Angstroms")

print("Vo+Deltao=",Vo, "eV")
print("m=",m,"eV")
print("hbar=",hbar,"eV*Angstrom")
print("e=",e,"Coulomb")
print("beta=",beta, "1/eV")

A=4*np.pi*np.sqrt(4*m)*dx/h #constante definida

area=346*375 #microns squared

conver=1E8#angstroms to microns 

eta=kb*T*np.log(3) #groundstate taken at 1/2 of Bose-Einsteins distribution

mu=eta#chemical potential, no particles are added


EV=np.linspace(0,1E-3,300)#(-50E15,10E15,1000)*2*e#applied potential in eV

V=EV/(2*e) #Applied potential in V


#Definiendo constantes de la integral usada

def C1():
    
    fact1=16.0*np.pi*m*e/(pow(h,3)*pow(beta,2)) #cambio hecho en beta
    
    fact2=np.exp(-A*np.sqrt(Vo))

    return fact1*fact2*hbar*conver/hs


def C2():

    return A/(2*np.sqrt(Vo)*pow(beta,1))# #cambio hecho en beta

def C3(EV):

    return np.exp((mu-EV)*beta)

def C4(EV):
    
    return np.exp(mu*beta)


#Definiendo la funcion a integrar

def func(x,EV):

    f1=np.exp(C2()*(x+eta*pow(beta,1)))# cambio hecho en beta
    f2=np.log( (1-C3(EV)*np.exp(-x))/ (1-C4(EV)*np.exp(-x))   )

    return C1()*f1*f2


J1=np.zeros(len(V))
J2=np.zeros(len(V))
u=np.ones(len(EV))*eta

print("C1=", C1())
print("C2=", C2())
print("C3=", C3(EV[0]))
print("C4=", C4(EV[0]))
print("A=",A )


for i in range(len(EV)):

    funcion=lambda x: func(x,EV[i])
    
    #resultados=integrate.quad(funcion,0,np.infty)
    
    #J1[i]=area*resultados[0]
    
    J2[i]=area*C1()*np.exp(beta*(C2()*eta+mu))*(1-np.exp(-beta*EV[i]))/(1-C2())



v=np.linspace(min(J2),max(J2),len(EV))

'''
    print("C1=", C1(a,m,h,beta,Vo,g))
    print("C2=", C2(a,m,Vo,beta))
    print("C3[",i,"]=", C3(mu,EV[i]))
    print("C4[",i,"]", C4(mu,EV[i]))


plt.plot(EV,J1)
plt.xlabel("$ V [Volts]  $")
plt.ylabel("$ J1 $")
plt.show()
plt.close()

'''

plt.plot(EV,J2)

plt.plot(u,v,"k--",linewidth=0.5)

plt.xlabel("$ V [Volts]  $",size=15)
plt.ylabel("$ I [Amp]$ (approx) ",size=15)
plt.title("$I-V\ curve$ ",size=15)
plt.savefig("IV_approx.png")
plt.show()
plt.close()















