import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad



alphaI=0.28125
alphaN=1.38
m=0
Io=0.904
B=0.51/2

mu=1 #a.u.
hbar=1 #a.u.
Z=2

Factor=24.18884 #time[as]/a.u.



def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)

def potential(r):
    
    Z=2.0
    ro=1.0/(2*Z)
    return -(Z-1)/r  + F*r*np.cos(theta) -(np.cos(theta)*alphaI*F/(r)**2)*np.exp(-3/r) +I(F)


def potential_schro(r):
    
    Z=2.0
    ro=1.0/(2*Z)
    return -(Z-1)/r -(1/r)*(1+r/(2*ro))*np.exp(-r/ro)  + F*r*np.cos(theta) -(np.cos(theta)*alphaI*F/(r)**2)*np.exp(-3/r) +I(F)

def kappa_C(r):

    return np.sqrt(2*mu*potential_schro(r))/hbar

def kappa(r):
    
    return np.sqrt(2*mu*potential(r))/hbar

def k():
    return np.sqrt(2*mu*I(F)/4)/hbar


r1=np.linspace(0.5,100,1000)

f1=[0.04,0.06,0.08,0.10,0.11]

theta=np.pi/3
THETA=[0,30,60,90,120,150,180]
for i in range(len(THETA)):
    F=0.11
    theta=THETA[i]*np.pi/180
    ones=np.ones(len(r1))*-I(F)
    func=potential_schro(r1)-I(F)
    plt.plot(r1,func,label=str(THETA[i]))
    plt.legend()
    plt.plot(r1,ones)

plt.ylim(-1,0.55)
plt.xlim(0,50)
plt.xlabel("r")
plt.ylabel("V(r)")
plt.title("Potential with different angles")
plt.savefig("potential_with_theta.eps")
plt.show()




#-----------------------------FINDING THE TURNING POINTS-----------------------------


f=np.linspace(0.04,0.11,8)
theta=4*np.pi/6
THETA=[1]
Colors=["k","b","g","c"]

#-----------------------------INTEGRATING-----------------------------



def inte_num(x):
    
    func=lambda y: kappa(y)
    return (np.exp(-2*quad(func,T1,x)[0]))

def inte_exp(T1,T2):
    
    func=lambda y: kappa(y)
    return np.exp(-2*quad(func,T1,T2)[0])


def inte_num_C(x):

    func=lambda y: kappa_C(y)
    return (np.exp(-2*quad(func,T1_C,x)[0]))

def inte_exp_C(T1_C,T2_C):
    
    func=lambda y: kappa_C(y)
    return np.exp(-2*quad(func,T1_C,T2_C)[0])


for j in range(len(THETA)):
    theta=THETA[j]*np.pi
    Time=[]
    Time_C=[]
    W=[]
    W_C=[]
    Turning_C=[]#Turning points of corrected function
    Turning=[]#Turning points of uncorrected function



    for i in f:
        F=i
        Turning_C.append([F,brentq(potential_schro,0.1,4),brentq(potential_schro,5,1000)])
    
        Turning.append([F,brentq(potential,0.1,4),brentq(potential,5,1000)])



#DWELL TIME___________________________________________________________


    for i in range(len(f)):
    
        F=f[i]
    
        T1=Turning[i][1]
        T2=Turning[i][2]
        W.append((T2-T1))
    
        T1_C=Turning_C[i][1]
        T2_C=Turning_C[i][2]
        W_C.append((T2_C-T1_C))
    

        func=lambda x: (1/kappa(x))*inte_num(x)
        func_C=lambda x: (1/kappa_C(x))*inte_num_C(x)
    
    
        exponencial=1#inte_exp(T1,T2)
        exponencial_C=1#inte_exp_C(T1_C,T2_C)
    
        Time.append(Factor*quad(func,T1,T2)[0]/exponencial)
        Time_C.append(Factor*quad(func_C,T1_C,T2_C)[0]/exponencial_C)

    print(j)




    plt.plot(f,Time,Colors[j],label="uncorrected, ")
    plt.plot(f,Time_C,Colors[j]+"--",label="corrected, ")

plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
#plt.ylim(0,200)
plt.title("Dwell time (in R), theta="+str(theta*180/np.pi))
plt.legend()
plt.savefig("dwell_time_R_"+str(theta*180/np.pi)+".eps")
plt.show()
plt.close()





















