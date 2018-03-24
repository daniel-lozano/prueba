import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate




dx=20.7989  # Angstrom
e=1.6E-19 # Coulomb
m=0.51E6   # eV
Do=30.5E-4 # eV
hbar=1973.27# eV*Angstrom =[hbar*c]
h=2*np.pi*hbar # eV*Angstrom
hs=6.5E-16 # eV s
kb=8.6E-5# eV/Kelvin/e
Tc=9.25 #Kelvin
g=(1.0/(92*2000))**2#(m/M)**2



print("dx=",dx, "Angstroms")
print("Deltao=",Do, "eV")
print("m=",m,"eV")
print("hbar=",hbar,"eV*Angstrom")
print("e=",e,"Coulomb")

A=4*np.pi*np.sqrt(4*m)*dx/h #constante definida

area=346*375 #microns squared

conver=1E8#angstroms to microns


EV=[0.0005,0.001,0.003,0.005]#(-50E15,10E15,1000)*2*e#applied potential in eV

#V=EV*(1/(2*e)) #Applied potential in V

T=np.linspace(0.1,1.0,10000)
J2=np.zeros(len(T))


eta=kb*T[0]*np.log(3) #groundstate taken at 1/2 of Bose-Einsteins distribution
mu=eta#chemical potential, no particles are added




#Definiendo constantes de la integral usada

def C1(Vo,beta):
    
    fact1=16.0*np.pi*m*e/(pow(h,3)*pow(beta,2)) #cambio hecho en beta
    
    fact2=np.exp(-A*np.sqrt(Vo))
    
    return fact1*fact2*hbar*conver/hs


def C2(Vo,beta):
    
    return A/(2*np.sqrt(Vo)*pow(beta,1))# #cambio hecho en beta

def C3(EV,Vo,beta):
    
    return np.exp((mu-EV)*beta)

def C4(EV,Vo,beta):
    
    return np.exp(mu*beta)


#Definiendo la funcion a integrar

def func(x,EV,beta):
    
    f1=np.exp(C2()*(x+eta*pow(beta,1)))# cambio hecho en beta
    f2=np.log( (1-C3(EV)*np.exp(-x))/ (1-C4(EV)*np.exp(-x))   )
    
    return C1()*f1*f2

#Definiendo funcion del gap

def Gap(T):
    
    return 1.8+3.07*0.5*Do*np.sqrt(1-T/Tc)

def DeT(T):
    
    if(T>0.6*Tc):
        
        return 3.07*kb*Tc*np.sqrt(1-T/Tc)
    
    if(T<=0.6*Tc):
        
        return factor*Do*(1-(2*np.pi*kb*T*np.exp(-Do/(kb*T)))/Do)
def Ns(T):
    
    mc2=2*0.51E6 #eV
    lam=47.0E-3 #nm
    TD=276.0 #k
    N0=19.87 # 1/eV
    g=1/(N0*np.log(1.13*TD/9.25))    #eV
    e=1.0#*1.6E-19
    
    ns=(mc2/(8*e*np.pi*pow(lam,2)))*(1-pow(T,4)*np.exp(4.0/(g*N0))/pow(1.13*TD,4))
    
    return ns

factor=0

a=3.07*kb*Tc*np.sqrt(1-0.6) #T>0.6Tc
    
b=Do*(1-(2*np.pi*kb*0.6*np.exp(-Do/(kb*0.6)))/Do)#T<=0.6Tc


factor=a/b
print(factor)

for k in range(len(EV)):

    for j in range(len(T)):
    
        beta=1.0/(kb*T[j]*Tc)
    
        eta=kb*T[j]*np.log(3)  #groundstate taken at 1/2 of Bose-Einsteins distribution
        
        mu=eta              #chemical potential, no particles are added
    
        l="$ T="+ str(T[j])+" K $"
    
        Vo=Gap(T[j]*Tc)
    
        DeltaT=DeT(T[j]*Tc)
        
        
        N=Ns(T[j]*Tc) #+ 0.5
        
        
        J2[j]=N*area*C1(Vo,beta)*np.exp(beta*(C2(Vo,beta)*eta+mu))*(1-np.exp(-beta*EV[k]))/(1-C2(Vo,beta))

    plt.plot(T,J2,label="$ V= $"+str(EV[k]))
    print(J2[-1])

plt.xlabel("$ T/T_c   $",size=15)
plt.ylabel("$ I [Amp]$ (approx) ",size=15)
plt.title("$I(T/T_c)\ curve$ ",size=15)
plt.legend(loc=2)

plt.savefig("I(T).png")
plt.show()
plt.close()



'''
T=np.linspace(0.8,1.0)
J2=np.zeros(len(T))

for k in range(len(EV)):
    
    for j in range(len(T)):
        
        beta=1.0/(kb*T[j]*Tc)
        
        eta=kb*T[j]*np.log(3)  #groundstate taken at 1/2 of Bose-Einsteins distribution
        
        mu=eta              #chemical potential, no particles are added
        
        l="$ T="+ str(T[j])+" K $"
        
        Vo=Gap(T[j]*Tc)
        
        DeltaT=DeT(T[j]*Tc)
        N=DeltaT/(2*g) +0.5#
        
        J2[j]=N*area*C1(Vo,beta)*np.exp(beta*(C2(Vo,beta)*eta+mu))*(1-np.exp(-beta*EV[k]))/(1-C2(Vo,beta))
    
    plt.plot(T,J2,label="$ V= $"+str(EV[k]))


plt.xlabel("$ T/T_c   $",size=15)
plt.ylabel("$ I [Amp]$ (approx) ",size=15)
plt.title("$I(T/T_c)\ curve$ ",size=15)
plt.legend(loc=3)

#plt.savefig("IV_approx_DT_lowT.png")
plt.show()
plt.close()



'''
















