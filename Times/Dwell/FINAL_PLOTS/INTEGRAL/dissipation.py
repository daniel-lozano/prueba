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
gamma=float(input("give a value of gamma= "))#float(argv[1])
print("gamma_D=",gamma)


f=np.linspace(0.04,0.11,7)

Factor=24.18884 #time[as]/a.u.

#------------------------------------------Energy factor-------------------------------------------
#D=0 #Dissipation

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2.)



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


#------------------------------------------Wave number without dissipation--------------------------------


def kappa_C(n):

    return np.sqrt(2*mu*abs(potential_schro(n)))/hbar

def kappa(n):
    
    return np.sqrt(2*mu*abs(potential(n)))/hbar

def k():
    
    return np.sqrt(2*mu*I(F)/4.)/hbar


#------------------------------------------Energy loss--------------------------------

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


def DE_func_num(T1,x):
    
    if(x>T1):
        func= lambda n: gamma*kappa(n)*hbar/mu
        return quad(func,T1,x)[0]
    if(x<=T1):
        return 0

def DE_func_num_C(T1_C,x):
    
    if(x>T1_C):
        func= lambda n: gamma*kappa_C(n)*hbar/mu
        return quad(func,T1_C,x)[0]
    if(x<=T1_C):
        return 0



#-------------------Dissipation potential--------------------------------


def disip_potential(n):
    
    return potential(n)+ DE_func_num(T1,n)

def disip_potential_C(n):
    
    return potential_schro(n) + DE_func_num_C(T1_C,n)

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

x=np.linspace(0.1,100)

Turning_C=[]#Turning points of corrected function
Turning=[]#Turning points of uncorrected function

print("F","T1","T2","T1_C","T2_C")

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
    
    print(F,Turning[-1][1],Turning[-1][2],Turning_C[-1][1],Turning_C[-1][2])


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




DE=np.zeros(len(f))
TE=np.zeros(len(f))

#----------Total Dissipative Energy---------------------------------

x_plot=np.linspace(0.1,100,1000)

for i in [0,1]:#range(len(f)):
    
    F=f[i]
    TE=np.ones(len(x_plot))*I(F)/4.0
    DE=DE_func(Turning_C[i][1],x_plot,Turning_C[i][2])
    plt.plot(x_plot,DE,label=str(F))
    plt.plot(x_plot,TE,"--",label="TE,F="+str(F))



plt.title("$ \mathrm{Dissipation\ in }\ \eta  $")
plt.xlabel("$ \eta\  \mathrm{(a.u.)} $")
plt.ylabel("$ \Delta E(\eta)\  \mathrm{(a.u.)} $")
plt.legend()
plt.savefig("Energy_lost.png")
plt.show()
plt.close()

#DWELL TIME___________________________________________________________



Time=[]
Time_C=[]
Time_ave=np.zeros(len(f))
Time_ave_C=np.zeros(len(f))
W=[]
W_C=[]
FILE=open("time_disip.txt","w")
FILE.write("#F Time Time_C \n")


for i in range(len(f)):
    
    F=f[i]
    
    T1=Turning[i][1]
    T2=Turning[i][2]
    W.append((T2-T1))
    
    T1_C=Turning_C[i][1]
    T2_C=Turning_C[i][2]
    W_C.append((T2_C-T1_C))
    
    #print(T1,T2)
    
    func=lambda x: (1.0/disip_kappa(x))*inte_num(x)
    func_C=lambda x: (1/disip_kappa_C(x))*inte_num_C(x)
    
    average_dt=Factor*quad(func,T1,T2)[0]
    average_dt_C=Factor*quad(func_C,T1_C,T2_C)[0]
    
    
    FILE.write(str(F)+" "+str(average_dt)+"\n")
    print(F,average_dt)
    Time_ave[i]=average_dt
    Time_ave_C[i]=average_dt_C
    
    exponencial=inte_exp(T1,T2)
    exponencial_C=inte_exp_C(T1_C,T2_C)
    Time.append(average_dt/exponencial)
    Time_C.append(average_dt_C/exponencial_C)
    print(F,exponencial,exponencial_C)
    
    
    if(i%1==0):
        print(i)

FILE.close()

print("Terminando la simulacion")

plt.plot(f,Time_ave,"k",label="uncorrected")
plt.plot(f,Time_ave_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Average dwell time, $ \\gamma= $"+str(gamma))
plt.legend()
plt.savefig("average_disip.png")
plt.show()
plt.close()


'''
plt.plot(f,Time,"k",label="uncorrected")
plt.plot(f,Time_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Transmission dwell time")
plt.legend()
plt.savefig("transmission_disip.png")
plt.show()
plt.close()

'''


















