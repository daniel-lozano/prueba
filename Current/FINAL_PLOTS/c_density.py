import numpy as np
import matplotlib.pyplot as plt
#from scipy import integrate




dx=20.7989  # Angstrom
e=1.6E-19 # Coulomb
m=0.51E6   # eV

Do=30.5E-4 # eV

hbar=1973.27# eV*Angstrom =[hbar*c]
h=2*np.pi*hbar # eV*Angstrom
hs=6.5E-16 # eV s

kb=8.6E-5# eV/Kelvin/e

Tc=9.25 #Kelvin
T=np.linspace(1,9,7)#9.2#4 Kelvin
t=T/Tc #Temperatura adimensional

g=(1.0/(92*2000))**2#(m/M)**2



print("dx=",dx, "Angstroms")
print("Delta_o=",Do, "eV")
print("m=",m,"eV")
print("hbar=",hbar,"eV*Angstrom")
print("e=",e,"Coulomb")
print("g=",g,"Coupling constant")

A=4*np.pi*np.sqrt(4*m)*dx/h #constante definida

area=1.0#346*375 #microns squared

conver=1E8#angstroms to microns 


EV=np.linspace(-0,0.5E-3,300)#/(2*Do)#(-50E15,10E15,1000)*2*e#applied potential in eV
print(EV)
V=EV/(2*e) #Applied potential in V



    
eta=0.0#kb*T[0]*np.log(3) #groundstate taken at 1/2 of Bose-Einsteins distribution
mu=0.0#chemical potential, no particles are added




#Definiendo constantes de la integral usada

def C1(Vo,beta):
    
    fact1=8.0*np.pi*pow(m,1)*e/(pow(h,3)*pow(beta,2)) #cambio hecho en beta
    
    fact2=np.exp(-A*np.sqrt(Vo))

    return fact1*fact2*hbar*conver/hs


def C2(Vo,beta):

    return A/(2*np.sqrt(Vo)*pow(beta,1))# #cambio hecho en beta

def C3(EV,Vo,beta):

    return np.exp((mu-EV*(2*Do))*beta)

def C4(EV,Vo,beta):
    
    return np.exp(mu*beta)


#Definiendo la funcion a integrar------------------------------------------

def func(x,EV,beta):

    f1=np.exp(C2()*(x+eta*pow(beta,1)))# cambio hecho en beta
    f2=np.log( (1-C3(EV)*np.exp(-x))/ (1-C4(EV)*np.exp(-x))   )

    return C1()*f1*f2

#Definiendo funcion del gap------------------------------------------------

def Gap(T):

    return 1.8+3.07*0.5*Do*np.sqrt(1-T/Tc)

def DeT(T):

    if(T>0.6*Tc):
        return 3.07*kb*Tc*np.sqrt(1-T/Tc)
    if(T<0.6*Tc):
        return Do*(1-(2*np.pi*kb*T*np.exp(-Do/(kb*T)))/Do)
#Definiendo funcion densidad de particulas----------------------------------

def Ns(T):
    
    mc2=2*0.51E6 #eV
    lam=47.0E-3 #nm
    TD=276.0 #k
    N0=19.87 # 1/eV
    g=1/(N0*np.log(1.13*TD/9.25))    #eV
    e=1.0#*1.6E-19

    ns=(mc2/(8.0*e*np.pi*pow(lam,2)))*(1-pow(T,4)*np.exp(4.0/(g*N0))/pow(1.13*TD,4))
    return ns


J1=np.zeros(len(V))
J2=np.zeros(len(V))
beta=1.0/(kb*T[0])

print("C1=", C1(Do,beta))
print("C2=", C2(Do,beta))
print("C3=", C3(EV[0],Do,beta))
print("C4=", C4(EV[0],Do,beta))
print("A=",A )

plt.figure(1,figsize=(15,5))
plt.subplot(121)
for j in range(len(T)):
    
    beta=1.0/(kb*T[j])

    eta=0#kb*T[j]*np.log(3) #groundstate taken at 1/2 of Bose-Einsteins distribution
    mu=eta#chemical potential, no particles are added
    
    u=np.ones(len(EV))*eta

    l="$ T/T_c="+ str(round(t[j],3))+"$"
    
    Vo=Gap(T[j])
    
    DeltaT=DeT(T[j])#3.07*kb*Tc*np.sqrt(1-T[j]/Tc)
    
    N=Ns(T[j]) #+0.5#
    print("N=",N)
    
    for i in range(len(EV)):
        #funcion=lambda x: func(x,EV[i])
    
        #resultados=integrate.quad(funcion,0,np.infty)
    
        #J1[i]=area*resultados[0]
        
        units=1#(2*DeT(T[j]))
        J2[i]=N*area*C1(Vo,beta)*np.exp(beta*(C2(Vo,beta)*eta+mu))*(1-np.exp(-beta*EV[i]*units))/(1-C2(Vo,beta))

    plt.plot(EV,J2,label=l)

#
#    v=np.linspace(min(J2),max(J2),len(EV))
#    plt.plot(u,v,"k--",linewidth=0.5)

caption="$ mu=k_b T ln(3)  $"
if(mu==0):
    caption="mu=0"
    
plt.xlabel("$ V  $",size=15)
plt.ylabel("$ J $ (approx) ",size=15)
plt.title("$ \mathrm{t}\ $ " + caption,size=15)
plt.legend(loc=2)


#plt.savefig("JV_approx_DT_lowT_mu="+str(round(mu,4))+".png")
#plt.show()
#plt.close()

plt.subplot(122)
T=np.linspace(9,9.249,5)
t=T/Tc
for j in range(len(T)):
    beta=1.0/(kb*T[j])
    
    eta=0#kb*T[j]*np.log(3) #groundstate taken at 1/2 of Bose-Einsteins distribution
    mu=eta#chemical potential, no particles are added
    
    u=np.ones(len(EV))*eta
    
    l="$ T/T_c="+ str(round(t[j],3))+"$"
    
    Vo=Gap(T[j])
    
    DeltaT=3.07*kb*Tc*np.sqrt(1-T[j]/Tc)
    
    N=Ns(T[j]) #+0.5#
    print("N=",N)
    
    for i in range(len(EV)):
        #funcion=lambda x: func(x,EV[i])
        
        #resultados=integrate.quad(funcion,0,np.infty)
        
        #J1[i]=area*resultados[0]
        units=1#(2*DeT(T[j]))
        J2[i]=N*area*C1(Vo,beta)*np.exp(beta*(C2(Vo,beta)*eta+mu))*(1-np.exp(-beta*EV[i]*units))/(1-C2(Vo,beta))

    plt.plot(EV,J2,label=l)


#v=np.linspace(min(J2),max(J2),len(EV))
#   plt.plot(u,v,"k--",linewidth=0.5)


caption="$ mu=k_b T ln(3)  $"
if(mu==0):
    caption="mu=0"
    
plt.xlabel("$ V  $",size=15)
plt.ylabel("$ J $ (approx) ",size=15)
plt.title("$ \mathrm{High\ t}\ $ " + caption,size=15)

plt.legend(loc=2)

plt.savefig("JV_approx_DT_mu="+str(round(mu,4))+".png")
plt.show()
plt.close()








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


















