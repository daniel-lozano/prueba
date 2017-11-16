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

Factor=24.18884

#F=0


def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)

def potential(n):
    
    
    ro=1/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(2*I(F))/2.0
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2)
    
    t5= 0
    
    return -t1 - t2 + t3 + t4 - t5 +I(F)/4.0


def potential_schro(n):
    
    
    ro=1/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(2*I(F))/2.0
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2)
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5 +I(F)/4.0

def kappa_C(n):

    return np.sqrt(2*mu*potential_schro(n))/hbar

def kappa(n):
    
    return np.sqrt(2*mu*potential(n))/hbar

'''
F=0.01
n=np.linspace(0.1,60,1000)
plt.plot(n,potential_schro(n))
plt.ylim(-0.5,0.5)
plt.show()
'''

'''
-----------------------------FINDING THE TURNING POINTS-----------------------------
'''

#f=np.linspace(0.1,0.8,15)
f=np.linspace(0.04,0.11,15)

Turning_C=[]
Turning=[]
FILE=open("turning_points.txt","w")
FILE.write("#F T1C T2C T1 T2 \n")
for i in f:
    F=i
    Turning_C.append([F,brentq(potential_schro,0.1,2),brentq(potential_schro,2,50)])
    
    Turning.append([F,brentq(potential,0.1,2),brentq(potential,2,50)])
    
    FILE.write(str(F)+" "+str(Turning_C[-1][1])+" "+str(Turning_C[-1][2])+" "+ str(Turning[-1][1])+" "+str(Turning[-1][2])+"\n")
    
    #print(F,Turning_C[-1][1],Turning_C[-1][2],Turning[-1][1],Turning[-1][2])
    
    #print(potential_schro(Turning_C[-1][1]),potential_schro(Turning_C[-1][2]),potential(Turning[-1][1]),potential(Turning[-1][2]),"\n")

FILE.close()
#print(Turning)

'''
-----------------------------INTEGRATING-----------------------------
'''


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



Time=[]
Time_C=[]
W=[]
W_C=[]


'''
DWELL TIME___________________________________________________________
'''



for i in range(len(f)):
    
    F=f[i]
    
    T1=Turning[i][1]
    T2=Turning[i][2]
    W.append(T2-T1)
    
    T1_C=Turning_C[i][1]
    T2_C=Turning_C[i][2]
    W_C.append(T2_C-T1_C)
    
    #print(T1,T2)
    
    func=lambda x: (1/kappa(x))*inte_num(x)
    func_C=lambda x: (1/kappa_C(x))*inte_num_C(x)
    
    
    exponencial=1#inte_exp(T1,T2)
    exponencial_C=1#inte_exp_C(T1_C,T2_C)
    
    Time.append(Factor*quad(func,T1,T2)[0]/exponencial)
    Time_C.append(Factor*quad(func_C,T1_C,T2_C)[0]/exponencial_C)



plt.plot(f,Time,"k",label="uncorrected")
plt.plot(f,Time_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Dwell time")
plt.legend()
plt.savefig("dwell_time.png")
plt.show()
plt.close()

'''
TRAVERSAL TIME___________________________________________________________
'''


T=[]
T_C=[]


for i in range(len(f)):
    
    F=f[i]
    
    T1=Turning[i][1]
    T2=Turning[i][2]
    
    T1_C=Turning_C[i][1]
    T2_C=Turning_C[i][2]
    
    #print(T1,T2)
    
    func=lambda x: (1/kappa(x))
    func_C=lambda x: (1/kappa_C(x))
    
    T.append(Factor*quad(func,T1,T2)[0])
    T_C.append(Factor*quad(func_C,T1_C,T2_C)[0])



plt.plot(f,T,"k",label="uncorrected")
plt.plot(f,T_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Traversal time")
plt.legend()
plt.savefig("traversal_time.png")
plt.show()
plt.close()


W_alex=[8.73,10.19,11.96,13.55,15.76,17.22,19.08,20.85]
T_alex=[34.30,39.83,46.29,51.85,60.15,65.68,73.04,79.50]


#plt.plot(W,T,label="Unconrrected")
#plt.plot(W_C,T_C,label="Corrected")
plt.plot(W_alex,T_alex,label="Alex")
plt.xlabel("Barrier width")
plt.ylabel("Time[as]")
plt.legend()
plt.show()
plt.close()

















