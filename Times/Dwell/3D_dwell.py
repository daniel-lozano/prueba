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

def k():
    return np.sqrt(2*mu*I(F)/4)/hbar



'''
-----------------------------FINDING THE TURNING POINTS-----------------------------
'''


f=np.linspace(0.04,0.11,10)
theta=np.linspace(0,np.pi/3,10)
zeros=np.zeros(len(f))

F_M,Theta_M=np.meshgrid(f,theta)

T_C1,T_1=np.meshgrid(zeros,zeros)
T_C2,T_2=np.meshgrid(zeros,zeros)


for i in range(len(f)):
    for j in range(len(theta)):
    
        F=f[i]*np.cos(theta[j])
    
        T_C1[i][j]=brentq(potential_schro,0.1,4)
        T_C2[i][j]=brentq(potential_schro,10,100)

        T_1[i][j]=brentq(potential,0.1,4)
        T_2[i][j]=brentq(potential,10,100)

    


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

T_3D,T_C3D=np.meshgrid(zeros,zeros)


for i in range(len(f)):
    
    for j in range(len(theta)):
        
        F=F_M[i][j]*np.cos(Theta_M[i][j])
    
        T1=T_1[i][j]
        T2=T_2[i][j]
        
    
        T1_C=T_C1[i][1]
        T2_C=T_C2[i][2]
        
        func=lambda x: (1/kappa(x))*inte_num(x)
        func_C=lambda x: (1/kappa_C(x))*inte_num_C(x)
    
    
        exponencial=1#inte_exp(T1,T2)
        exponencial_C=1#inte_exp_C(T1_C,T2_C)
    
        T_3D[i][j]=Factor*quad(func,T1,T2)[0]/exponencial
        T_C3D[i][j]=Factor*quad(func_C,T1_C,T2_C)[0]/exponencial_C

    print(i)



"""
plt.plot(f,Time,"k",label="uncorrected")
plt.plot(f,Time_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Dwell time")
plt.legend()
plt.savefig("dwell_time.png")
plt.show()
plt.close()
"""



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot a basic wireframe.
ax.plot_wireframe(F_M, Theta_M, T_C3D, rstride=10, cstride=10)

plt.show()
