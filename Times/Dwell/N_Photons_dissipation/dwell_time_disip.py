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
print("omega=",omega)



Factor=24.18884 #time[as]/a.u.

#------------------------------------------Energy factor-------------------------------------------
D=0 #Dissipation

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)

def I_mod(F):
    
    return I(F)-D

#------------------------------------------Original potential factor-------------------------------------------

def potential(n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*pow(n,2))
    
    t4=  (alphaI*F/n**2) * np.exp(-3/n)
    
    
    return -t1 - t2 + t3 + t4  +I_mod(F)/4.0


def potential_schro(n):
    
    ro=1.0/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=   (alphaI*F/n**2) * np.exp(-3/n)
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    return -t1 - t2 + t3 + t4 - t5 +I_mod(F)/4.0

#------------------------------------------Wave number-------------------------------------------


def kappa_C(n):

    return np.sqrt(2*mu*abs(potential_schro(n)))/hbar

def kappa(n):
    
    return np.sqrt(2*mu*abs(potential(n)))/hbar

def k():
    return np.sqrt(2*mu*I(F)/4)/hbar

def DE_func(T1,T2):

    func= lambda n: gamma*kappa(n)*hbar/mu
    return quad(func,T1,T2)[0]



#-------------------Dissipation--------------------------------


def disip_potential(n):
    
    return potential(n)-n*gamma*(kappa(n)*hbar/(mu*8))

def disip_potential_C(n):
    
    return potential_schro(n)-n*gamma*(kappa_C(n)*hbar/(mu*8))

def disip_kappa_C(n):
    
    return np.sqrt(2*mu*abs(disip_potential_C(n)))/hbar

def disip_kappa(n):
    
    return np.sqrt(2*mu*abs(disip_potential(n)))/hbar




#encontrar el punto de retorno
def find(func,x):
    for i in range(len(func)):
        if(func[i]>0):
            return x[i]


'''
-----------------------------FINDING THE TURNING POINTS-----------------------------
'''

#f=np.linspace(0.1,0.8,15)
f=np.linspace(0.04,0.11,10)
x=np.linspace(0.1,20)

Turning_C=[]#Turning points of corrected function
Turning=[]#Turning points of uncorrected function
FILE=open("turning_points.txt","w")
FILE.write("#F T1C T2C T1 T2 \n")


for i in f:
    F=i
    
    keldish=omega*np.sqrt(2*(I(F)))/F
    print("F=",round(F,3),"gamma_k=",round(keldish,3))
    
    
    RETURN1=find(potential_schro(x),x)#2
    RETURN2=find(potential(x),x)#2
    
    Turning_C.append([F,brentq(potential_schro,0.1,RETURN1),brentq(potential_schro,RETURN1,100)])
    
    Turning.append([F,brentq(potential,0.1,RETURN2),brentq(potential,RETURN2,100)])
    
    FILE.write(str(F)+" "+str(Turning_C[-1][1])+" "+str(Turning_C[-1][2])+" "+ str(Turning[-1][1])+" "+str(Turning[-1][2])+"\n")
    
    #posibles errores grandes
    
    if(abs(potential_schro(Turning_C[-1][1])>1E-5)):
        print("Turning problem! in C1",str(F), potential_schro(Turning_C[-1][1]))
    
    if(abs(potential_schro(Turning_C[-1][2])>1E-5)):
        print("Turning problem! in C2",str(F),potential_schro(Turning_C[-1][2]))

    if(abs(potential(Turning[-1][1])>1E-5)):
        print("Turning problem! in 1",str(F),potential(Turning[-1][1]))

    if(abs(potential(Turning[-1][2])>1E-5)):
        print("Turning problem! in 2",str(F),potential(Turning[-1][2]))



FILE.close()
#print(Turning)


#-----------------------------INTEGRATING-----------------------------



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
DE=np.zeros(len(f))
TE=np.zeros(len(f))

#----------Total Dissipative Energy---------------------------------

#gamma=1E-2#E-1#4.5E-1
for j in range(4):
    
    for i in range(len(f)):
        
        gamma=(j+1.0)*3.0E-3
        F=f[i]
        DE[i]=DE_func(Turning_C[i][1],Turning_C[i][2])
        TE[i]=I(F)/4.0

#Finding turning points with dissipation------------------------------

Turning_C=[]#Turning points of corrected function
Turning=[]#Turning points of uncorrected function
cutoff=0.5
gamma=0.0005

print "gamma=",gamma
for i in range(len(f)):
    F=f[i]
    
    N=int(DE[i]/(hbar*omega))
    
    if(DE[i]/(hbar*omega)> N+cutoff):
        N+=1
    
    
    D=DE[i]#N*hbar*omega#
    
    print "F=",F, "N=", N, DE[i]/(hbar*omega) #campo y numero de fotones absorvidos
    
    
    RETURN1=find(potential_schro(x),x)#2
    RETURN2=find(potential(x),x)#2
    
    Turning_C.append([F,brentq(potential_schro,0.1,RETURN1),brentq(potential_schro,RETURN1,100)])
    
    Turning.append([F,brentq(potential,0.1,RETURN2),brentq(potential,RETURN2,100)])
    
    #posibles errores grandes
    
    if(abs(potential_schro(Turning_C[-1][1])>1E-5)):
        print("Turning problem! in C1",str(F), potential_schro(Turning_C[-1][1]))
    
    if(abs(potential_schro(Turning_C[-1][2])>1E-5)):
        print("Turning problem! in C2",str(F),potential_schro(Turning_C[-1][2]))
    
    if(abs(potential(Turning[-1][1])>1E-5)):
        print("Turning problem! in 1",str(F),potential(Turning[-1][1]))
    
    if(abs(potential(Turning[-1][2])>1E-5)):
        print("Turning problem! in 2",str(F),potential(Turning[-1][2]))




#DWELL TIME___________________________________________________________


#'''
for i in range(len(f)):
    
    F=f[i]
    
    N=int(DE[i]/(hbar*omega))
    
    if(DE[i]/(hbar*omega)> N+cutoff):
        N+=1
    
    D=DE[i]#N*hbar*omega#DE[i]
    

    T1=Turning[i][1]
    T2=Turning[i][2]
    W.append((T2-T1)/1.0)
    
    T1_C=Turning_C[i][1]
    T2_C=Turning_C[i][2]
    W_C.append((T2_C-T1_C)/1.0)
    
    #print(T1,T2)
    
    func=lambda x: (1/disip_kappa(x))*inte_num(x)
    func_C=lambda x: (1/disip_kappa_C(x))*inte_num_C(x)
    
    
    #exponencial=inte_exp(T1,T2)
    #exponencial_C=inte_exp_C(T1_C,T2_C)
    #print(F,exponencial,exponencial_C)
    
    Time_ave.append(Factor*quad(func,T1,T2)[0])
    Time_ave_C.append(Factor*quad(func_C,T1_C,T2_C)[0])
    #Time.append(Factor*quad(func,T1,T2)[0]/exponencial)
    #Time_C.append(Factor*quad(func_C,T1_C,T2_C)[0]/exponencial_C)
    
    if(i%1==0):
        print "Corrida=",i,"F=",F,"N=",N





plt.plot(f,Time_ave,"k-",label="uncorrected")
plt.plot(f,Time_ave_C,"r-",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Average dwell time, $ \gamma= $"+str(gamma))
plt.legend()
plt.savefig("gamma="+str(gamma)+"_average_disip.png")
#plt.show()
#plt.close()

'''#

plt.plot(f,Time,"k",label="uncorrected")
plt.plot(f,Time_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Transmission dwell time")
plt.legend()
plt.savefig("transmission_disip.png")
plt.show()
plt.close()

'''#


















