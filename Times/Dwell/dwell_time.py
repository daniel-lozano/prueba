import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad


alphaI=7.2
alphaN=11.08
m=0
Io=0.58
B=0.51/2

mu=0.0005484
hbar=1#0.197E9


F=0

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)



def potential_schro(n):
    
    Z=2
    
    ro=1/(2*Z)
    
    t1= (1-(1+m)*np.sqrt(I(F))/2.0)/(2*n)#(2-(1+m)*np.sqrt(I(F))/2.0)/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=  np.exp(-3/n) * (alphaI*F/n**2) +I(F)/4.0
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5

def kappa(n):

    return np.sqrt(2*mu*potential_schro(n))/hbar

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

f=np.linspace(0.04,0.12,15)

Turning=[]

for i in f:
    F=i
    Turning.append([F,brentq(potential_schro,0.1,4),brentq(potential_schro,4,30)])
    print(Turning[-1][1],Turning[-1][2])
    print(potential_schro(Turning[-1][1]),potential_schro(Turning[-1][2]),"\n")


#print(Turning)

'''
-----------------------------INTEGRATING-----------------------------
'''




Time=[]

def expo(n):
    
    func= lambda x: kappa(x)
    
    p1=0
    
    Res=[]
    for i in range(len(n)):
        print(quad(func,p1,n)[0])
        
        #Res.append(quad(funci,p1,n)[0])
    
#return Res#np.exp(-2*quad(funci,T1,n)[0])

print(expo([0,1,2,3,4]))

T1=Turning[1][1]
#n=np.linspace(0.1,20)
#plt.plot(n,expo(n))
#plt.show()


'''

for i in range(len(f)):
    
    F=f[i]
    
    T1=Turning[i][1]
    T2=Turning[i][2]
    
    print(T1,T2)
    
    func=lambda n: expo(n,T1)/kappa(n)
    
    Time.append(quad(func,T1,T2)[0])

print(Time)


'''












