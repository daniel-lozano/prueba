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
    
    
    return -t1 - t2 + t3 + t4  #+I(F)/4.0


def potential_schro(n):
    
    ro=1/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=   (alphaI*F/n**2) * np.exp(-3/n)
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5 +I(F)/4.0


def triangle(n):
    
    y=potential(n)
    
    m=max(y)
    
    maxIndexList = [index for index,value in enumerate(y) if value==max(y)]
    index=maxIndexList[0]

    
    no=n[index]
    ym=potential(no)
    
    y1=0
    y2=0
    x1=0
    x2=0
    
    function=np.zeros(len(n))
    
    up=4
    
    for i in range(len(n)):
        
       
        if(n[i]<up):
            
            
            x1=1.4
            x2=up
            y1=potential(x1)
            y2=potential(x2)
            
           
        
        
            
        if(n[i]>=up):
            x1=10
            x2=50
            y1=potential(x1)
            y2=potential(x2)
            
        M=(y2-y1)/(x2-x1)
        y0=y1-M*x1


        function[i]=M*n[i] +y0
            
            
            
    return function



'''
-----------------------------Ploting-----------------------------
'''
f=[0.02,0.04,0.06,0.08]
n=np.linspace(0.1,50)
marker=["b","r","g","c"]

for i in range(len(f)):
    
    F=f[i]
    V_ori=potential(n)
    V_tri=triangle(n)
    plt.plot(n,V_ori,marker[i],label=str(f[i]))
    plt.plot(n,V_tri,marker[i])



plt.ylim(-0.25,0)
plt.xlim(0,max(n))
plt.legend()
plt.show()
plt.close()











