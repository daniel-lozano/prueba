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

def potential(n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*pow(n,2))
    
    t4=  (alphaI*F/n**2) * np.exp(-3/n)
    
    
    return -t1 - t2 + t3 + t4  +I(F)/4.0


def potential_schro(n):
    
    ro=1/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=   (alphaI*F/n**2) * np.exp(-3/n)
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5 +I(F)/4.0

def kappa_C(n):
    
    return np.sqrt(2*mu*potential_schro(n))/hbar

def kappa(n):
    
    return np.sqrt(2*mu*potential(n))/hbar

def kf():
    return np.sqrt(2*mu*I(F)/4)/hbar
def kf_v(max):
    return np.sqrt(2*mu*abs(max))/hbar


'''
    -----------------------------FINDING THE TURNING POINTS-----------------------------
    '''

#f=np.linspace(0.1,0.8,15)
f=np.linspace(0.04,0.11,10)

Turning_C=[]#Turning points of corrected function
Turning=[]#Turning points of uncorrected function

for i in f:
    
    F=i
    
    Turning_C.append([F,brentq(potential_schro,0.1,2),brentq(potential_schro,10,50)])
    
    Turning.append([F,brentq(potential,0.1,2),brentq(potential,10,50)])
    

'''
    -----------------------------INTEGRATING-----------------------------
    '''

W=[]
W_C=[]
Max=[]
Max_C=[]
n=np.linspace(0.01,0.12,10)

for i in range(len(f)):
    
    F=f[i]
    
    T1=Turning[i][1]
    T2=Turning[i][2]
    W.append(T2-T1)
    Max.append(max(potential(n))+1)
    
    T1_C=Turning_C[i][1]
    T2_C=Turning_C[i][2]
    W_C.append(T2_C-T1_C)
    Max_C.append(max(potential_schro(n))+1)


    
Time_E=[]
Time_C=[]

for i in range(len(f)):
    
    F=f[i]
    a=W[i]
    k=kf()
    p=kf_v(Max[i])

    a=W[i]

    T_res=2*mu*(k*(p**2-k**2)*a+(k*(p**2+k**2)/(2*p))*np.sinh(2*p*a))/ (((p**2+k**2)*np.cosh(p*a))**2 -(p**2-k**2)**2 )
    Time_E.append(Factor*T_res)


plt.plot(f,Time_E)
plt.ylabel("time")
plt.xlabel("Field(a.u)")
plt.show()
plt.close()




















