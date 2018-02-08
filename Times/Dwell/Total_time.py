
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad



alphaI=0.28125
alphaN=1.38
m=0
Io=0.904
B=0.51/2
VALUE=2.4E-5
print("VALUE=",VALUE)
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
    if(r>VALUE):
        return -(Z-1)/r -(1/r)*(1+r/(2*ro))*np.exp(-r/ro)  + F*r*np.cos(theta) -(np.cos(theta)*alphaI*F/(r)**2)*np.exp(-3/r) +I(F)
    if(r<=VALUE):
        return -(Z)/VALUE +1/r-(1/r)*(1+r/(2*ro))*np.exp(-r/ro)  + F*r*np.cos(theta) -(np.cos(theta)*alphaI*F/(r)**2)*np.exp(-3/r) +I(F)

def kappa_C(r):
    
    return np.sqrt(2*mu*abs(potential_schro(r)))/hbar

def kappa(r):
    
    return np.sqrt(2*mu*abs(potential(r)))/hbar

def k():
    return np.sqrt(2*mu*I(F)/4)/hbar




#-----------------------------FINDING THE TURNING POINTS-----------------------------


f=np.linspace(0.04,0.11,8)
theta=95*np.pi/180


Turning_C=[]#Turning points of corrected function
Turning=[]#Turning points of uncorrected function

FILE=open("turning_points_R.txt","w")
FILE.write("#F T1C T2C T1 T2 \n")

for i in f:
    F=i
    Turning_C.append([F,brentq(potential_schro,0.1,4),brentq(potential_schro,5,50)])
    
    Turning.append([F,brentq(potential,0.1,4),brentq(potential,5,50)])
    
    FILE.write(str(F)+" "+str(Turning_C[-1][1])+" "+str(Turning_C[-1][2])+" "+ str(Turning[-1][1])+" "+str(Turning[-1][2])+"\n")


FILE.close()




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



Time=[]
Time_C=[]
W=[]
W_C=[]


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
    
    print(i)


plt.plot(f,Time,"k",label="uncorrected")
plt.plot(f,Time_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
#plt.ylim(0,200)
plt.title("Dwell time (in R), theta="+str(theta*180/np.pi))
plt.legend()
#plt.savefig("dwell_time_R_"+str(theta*180/np.pi)+".png")
plt.show()
plt.close()


#TRAVERSAL TIME___________________________________________________________



T=[]
T_C=[]


for i in range(len(f)):
    
    F=f[i]
    
    T1=0
    T2=Turning[i][1]
    
    T1_C=0
    T2_C=Turning_C[i][1]
    
    
    
    func=lambda x: (1/kappa(x))
    func_C=lambda x: (1/kappa_C(x))
    
    T.append(2*Factor*quad(func,T1,T2)[0])
    T_C.append(2*Factor*quad(func_C,T1_C,T2_C)[0])


plt.plot(f,T,"k",label="uncorrected")
plt.plot(f,T_C,"k--",label="corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.title("Traversal time (in R)")
plt.legend()
plt.savefig("traversal_time_R_region1.png")
plt.show()
plt.close()


W_alex=[8.73,10.19,11.96,13.55,15.76,17.22,19.08,20.85]
T_alex=[34.30,39.83,46.29,51.85,60.15,65.68,73.04,79.50]

plt.plot(W,f,label="Width with uncorrected")
plt.plot(W_C,f,label="Width with corrected")
plt.xlabel("Barrier width")
plt.ylabel("Field")
plt.legend()
plt.show()
plt.close()



plt.plot(W,Time,label="Unconrrected")
plt.plot(W_C,Time_C,label="Corrected")
plt.plot(W_alex,T_alex,label="Alex")
plt.xlabel("Barrier width")
plt.ylabel("Time[as]")
plt.legend()
plt.show()
plt.close()

#'''



















