import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad

from sys import argv



alphaI=0.28125
alphaN=1.38
m=0
Io=0.904#
B=0.51/2

mu=1.0 #a.u.
hbar=1.0 #a.u.

h=2.0*np.pi
c=137.0
lam=735.0/0.0529

Z=2
omega=h*c/lam #
number=0.0#float(argv[1])
gamma=float(argv[1])
print("gamma_D=",gamma)



Factor=24.18884 #time[as]/a.u.

#------------------------------------------Energy factor-------------------------------------------
#D=0 #Dissipation

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)



#------------------------------------------Original potential factor-------------------------------------------

def potential(n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*pow(n,2))
    
    t4=  (alphaI*F/n**2) * np.exp(-3/n)
    
    
    return -t1 - t2 + t3 + t4  +I(F)/4.0


def potential_schro(n):
    
    ro=1.0/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=   (alphaI*F/n**2) * np.exp(-3/n)
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    return -t1 - t2 + t3 + t4 - t5 +I(F)/4.0


#------------------------------------------Wave number withoud dissipation--------------------------------


def kappa_C(n):

    return np.sqrt(2*mu*abs(potential_schro(n)))/hbar

def kappa(n):
    
    return np.sqrt(2*mu*abs(potential(n)))/hbar

def k():
    return np.sqrt(2*mu*I(F)/4)/hbar

def DE_func(T1,x,T2):
    
    RES=np.zeros(len(x))
    
    for i in range(len(x)):
        func= lambda n: gamma*kappa(n)*hbar/mu
       
        if(x[i]<=T1):
           RES[i]=0
               
        if(x[i]>T1):
            RES[i]=quad(func,T1,x[i])[0]
        
        if(x[i]>T2):
            RES[i]=quad(func,T1,T2)[0]
        


    return RES

def DE_func_num(T1,x):
    
    
    if(x>T1):
        func= lambda n: gamma*kappa(n)*hbar/mu
        return quad(func,T1,x)[0]
    if(x<=T1):
        return 0


def DE_func_C(T1_C,x,T2_C):
    
    RES_C=np.zeros(len(x))
    
    
    
    for i in range(len(x)):
        func= lambda n: gamma*kappa_C(n)*hbar/mu
        
        if(x[i]<=T1_C):
            RES[i]=0
        
        if(x[i]>T1_C):
            RES_C[i]=quad(func,T1_C,x[i])[0]
        
        if( x[i]>T2_C):
            RES_C[i]=quad(func,T1_C,T2_C)[0]


    return RES_C

def DE_func_num_C(T1_C,x):
    
    if(x>T1_C):
        func= lambda n: gamma*kappa_C(n)*hbar/mu
        return quad(func,T1_C,x)[0]
    if(x<=T1_C):
        return 0



#-------------------Dissipation potential--------------------------------


def disip_potential(n):
    
    return potential(n)+DE_func_num(T1,n)

def disip_potential_C(n):
    
    return potential_schro(n)+ DE_func_num_C(T1_C,n)

#---------------- Dissipative potential--------------------------

def disip_kappa(n):
    
    return np.sqrt(2*mu*(disip_potential(n)))/hbar

def disip_kappa_C(n):
    
    return np.sqrt(2*mu*(disip_potential_C(n)))/hbar


#-----------------Funcion para encontrar el punto de retorno ------
def find(func,x):
    for i in range(len(func)):
        if(func[i]>0):
            return x[i]



#-----------------------------FINDING THE FIRST TURNING POINTS-----------------------------

#f=np.linspace(0.1,0.8,15)
f=np.linspace(0.04,0.11,7)
x=np.linspace(0.1,100)

Turning_C=[]#Turning points of corrected function
Turning=[]#Turning points of uncorrected function

for i in f:
    F=i
    
    RETURN1_C=find(potential_schro(x),x)#2
    RETURN1=find(potential(x),x)#2
    
    Turning_C.append([F,brentq(potential_schro,0.1,RETURN1_C),0])
    
    Turning.append([F,brentq(potential,0.1,RETURN1),0])
    
    #first turning point
    T1=Turning[-1][1]
    T1_C=Turning_C[-1][1]

    #finding second turning point
    Turning[-1][2]=brentq(disip_potential,RETURN1,100)
    Turning_C[-1][2]=brentq(disip_potential_C,RETURN1_C,100)
    
    print(F,Turning[-1][1],Turning[-1][2])

x_plot=np.linspace(0.5,x[-1],100)
'''
for i in range(len(f)):
    pos=i
    F=f[pos]
    

    potential_plot=potential(x_plot)-I(F)/4.0
    energy_plot=-I(F)/4.0-DE_func(Turning[pos][1],x_plot,Turning[pos][2])
    plt.plot(x_plot,potential_plot)
    plt.plot(x_plot,energy_plot)
    plt.title(str(F))
    plt.xlim(0,x_plot,Turning[pos][2]+1)

    plt.show()

    plt.close()
'''





print("Se encontraron los puntos de retorno")

#-----------------------------INTEGRATING-----------------------------
#'''


def inte_num(x):
    
    func=lambda y: disip_kappa(y)
    return (np.exp(-2*quad(func,T1,x)[0]))

def inte_exp(T1,T2):
    
    func=lambda y: disip_kappa(y)
    return np.exp(-2*quad(func,T1,T2)[0])


def inte_num_C(x):

    func=lambda y: disip_kappa_C(y)
    return (np.exp(-2*quad(func,T1_C,x)[0]))

def inte_exp_C(T1_C,T2_C):
    
    func=lambda y: disip_kappa_C(y)
    return np.exp(-2*quad(func,T1_C,T2_C)[0])



Time=[]
Time_C=[]
Time_ave=[]
Time_ave_C=[]
W=[]
W_C=[]


#----------Total Dissipative Energy---------------------------------

LAB=["k","b","r","g","c","y","m"]

for i in range(len(f)):
    pos=i
    F=f[pos]
    TE=np.ones(len(x_plot))*(-I(F)/4.0)-DE_func(Turning[pos][1],x_plot,Turning[pos][2])
    POT=potential(x_plot)-I(F)/4.0
    plt.plot(x_plot,POT,LAB[i],label="F="+str(F))
    plt.plot(x_plot,TE,LAB[i]+"--")
    plt.legend(loc=1)
    plt.xlim(0,80)
    plt.ylim(-0.5,0)
    plt.xlabel("$\eta\ \mathrm{(a.u.)} $")
    plt.ylabel("$ V(\eta),\ E(\eta)\ \mathrm{(a.u.)} $")

plt.savefig("true_system.png")
plt.show()
plt.close()


























