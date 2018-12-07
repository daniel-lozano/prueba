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

mu=1 #a.u.
hbar=1 #a.u.

h=2.0*np.pi
c=137.0
lam=735.0/0.0529

Z=2
omega=h*c/lam#
number=0
print("omega=",omega)

Factor=24.18884 #time[as]/a.u.

def I(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2)

def I_mod(F):
    
    return Io + ((alphaN-alphaI)*F**2)/(2) - omega*number

def potential(n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*pow(n,2))
    
    t4=  (alphaI*F/n**2) * np.exp(-3/n)
    
    
    return -t1 - t2 + t3 + t4  +I_mod(F)/4.0

# Expressions for the triangle potential, the derivative is define for the slope of the triangle in the second turning point. The potential will be define as a general expression and not by parts since the expressions used will be define in the domain which its valid.

def potential_derivative(n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n**2)
    
    t2= F/8
    
    t3= (m**2-1)/(4*pow(n,3))
    
    t4=  (alphaI*F/n**2) *( np.exp(-3/n)*(3-2*n)/n**4  )

    return +t1 - t2 - t3 + t4

# The dependence on the energy is necessary to determine the point to take the slope as well as the turning point, this is taken into account in the original expression as it has the term needed

def potential_triangle(n,n1,n2):

    m=potential_derivative(n2) #Slope of the triangle
    b=-m*n2 # cut with the axis
    
    return m*n+b



def potential_schro(n):
    
    ro=1.0/(2*Z)
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n)
    
    t2= n*F/8
    
    t3= (m**2-1)/(8*n**2)
    
    t4=   (alphaI*F/n**2) * np.exp(-3/n)
    
    t5= (1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    
    return -t1 - t2 + t3 + t4 - t5 +I_mod(F)/4.0

# Expressions for the triangle potential, the derivative is define for the slope of the triangle in the second turning point. The potential will be define as a general expression and not by parts since the expressions used will be define in the domain which its valid.

def potential_derivative_C(n):
    
    b2=(Z-1)-(1+m)*np.sqrt(I(F)/2)
    
    t1= b2/(2*n**2)
    
    t2= F/8
    
    t3= (m**2-1)/(4*pow(n,3))
    
    t4=  (alphaI*F/n**2) *( np.exp(-3/n)*(3-2*n)/n**4  )
    
    t5= -(1/n**2)*np.exp(-n/(2*ro)) - (1/(2*ro))*(1/n +1/(4*ro))*np.exp(-n/(2*ro))
    
    return +t1 - t2 - t3 + t4 - t5

# The dependence on the energy is necessary to determine the point to take the slope as well as the turning point, this is taken into account in the original expression as it has the term needed

def potential_triangle_C(n,n1,n2):
    
    m=potential_derivative_C(n2) #Slope of the triangle
    b=-m*n2 # cut with the axis
    
    return m*n+b



def kappa_C(n):

    return np.sqrt(2*mu*abs(potential_schro(n)))/hbar

def kappa(n):
    
    return np.sqrt(2*mu*abs(potential(n)))/hbar

def kappa_T(n,n1,n2):
    
    return np.sqrt(2*mu*abs(potential_triangle(n,n1,n2)))/hbar

def kappa_T_C(n,n1,n2):
    
    return np.sqrt(2*mu*abs(potential_triangle_C(n,n1,n2)))/hbar



def k():
    return np.sqrt(2*mu*I(F)/4)/hbar

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
    
    if(abs(potential_schro(Turning_C[-1][1])>1E-10)):
        print("Turning problem! in C1",str(F), potential_schro(Turning_C[-1][1]))
    
    if(abs(potential_schro(Turning_C[-1][2])>1E-10)):
        print("Turning problem! in C2",str(F),potential_schro(Turning_C[-1][2]))

    if(abs(potential(Turning[-1][1])>1E-10)):
        print("Turning problem! in 1",str(F),potential(Turning[-1][1]))

    if(abs(potential(Turning[-1][2])>1E-10)):
        print("Turning problem! in 2",str(F),potential(Turning[-1][2]))



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

# Integrals for the triangle potential

def inte_num_T(x,n1,n2):
    
    func=lambda y: kappa_T(y,n1,n2)
    return (np.exp(-2*quad(func,T1,x)[0]))

def inte_exp_T(T1,T2,n1,n2):
    
    func=lambda y: kappa_T(y,n1,n2)
    return np.exp(-2*quad(func,T1,T2)[0])


# Integrals for the corrected potential

def inte_num_C(x):

    func=lambda y: kappa_C(y)
    return (np.exp(-2*quad(func,T1_C,x)[0]))

def inte_exp_C(T1_C,T2_C):
    
    func=lambda y: kappa_C(y)
    return np.exp(-2*quad(func,T1_C,T2_C)[0])

# Integrals for the triangle potential

def inte_num_T_C(x,n1,n2):
    
    func=lambda y: kappa_T_C(y,n1,n2)
    return (np.exp(-2*quad(func,T1,x)[0]))

def inte_exp_T_C(T1,T2,n1,n2):
    
    func=lambda y: kappa_T_C(y,n1,n2)
    return np.exp(-2*quad(func,T1,T2)[0])






Time1=[]
Time_C1=[]
Time_T1=[]
Time_TC1=[]

Time2=[]
Time_C2=[]
Time_T2=[]
Time_TC2=[]



W=[]
W_C=[]


'''
DWELL TIME___________________________________________________________
'''



for i in range(len(f)):
    
    F=f[i]
    
    T1=Turning[i][1]
    T2=Turning[i][2]
    W.append((T2-T1)/1.0)
    
    T1_C=Turning_C[i][1]
    T2_C=Turning_C[i][2]
    W_C.append((T2_C-T1_C)/1.0)
    
    #print(T1,T2)
    
    func=lambda x: (1/kappa(x))*inte_num(x)
    func_T=lambda x: (1/kappa_T(x,T1,T2))*inte_num_T(x,T1,T2)
    func_C=lambda x: (1/kappa_C(x))*inte_num_C(x)
    func_T_C=lambda x: (1/kappa_T_C(x,T1,T2))*inte_num_T_C(x,T1_C,T2_C)

    
    
    exponencial=1#inte_exp(T1,T2)
    exponencial_T=1#inte_exp_T(T1,T2)
    exponencial_C=1#inte_exp_C(T1_C,T2_C)
    exponencial_T_C=1#inte_exp_T_C(T1_C,T2_C)

    
    Time2.append(Factor*quad(func,T1,T2)[0]/exponencial)
    Time_T2.append(Factor*quad(func_T,T1,T2)[0]/exponencial_T)
    Time_C2.append(Factor*quad(func_C,T1_C,T2_C)[0]/exponencial_C)
    Time_TC2.append(Factor*quad(func_T,T1_C,T2_C)[0]/exponencial_T_C)

    
    #turning points
    T1=0
    T2=Turning[i][1]
    T1_C=0
    T2_C=Turning_C[i][1]
    
    exponencial=1#inte_exp(T1,T2)
    exponencial_C=1#inte_exp_C(T1_C,T2_C)
    
    Time1.append(Factor*quad(func,T1,T2)[0]/exponencial)
    Time_C1.append(Factor*quad(func_C,T1_C,T2_C)[0]/exponencial_C)
    
    if(i%2==0):
        print(i)
#print("dif",abs(Time2[-1]-Time_C2[-1]))



plt.subplot(211)
plt.plot(f,Time1,"b",label="uncorrected Region 1")
plt.plot(f,Time_C1,"b--",label="corrected Region 1")
plt.ylabel("Time [as]")
plt.title("Dwell time ")
plt.legend()
plt.subplot(212)
plt.plot(f,Time2,"k",label="Uncorrected Region 2")
plt.plot(f,Time_C2,"k--",label="Corrected Region 2")
plt.plot(f,Time_T2,"r",label="triangle Uncorrected Region 2")
plt.plot(f,Time_TC2,"r--",label="triangle Corrected Region 2")

plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.legend()
plt.savefig("dwell_time.png")
plt.show()
plt.close()


'''
TRAVERSAL TIME___________________________________________________________
'''


T_1=[]
T_C1=[]
T_2=[]
T_C2=[]


for i in range(len(f)):
    
    F=f[i]
    
    T1=0#Turning[i][1]
    T2=Turning[i][1]#
    
    T1_C=0#Turning_C[i][1]
    T2_C=Turning_C[i][1]#
    
    #print(T1,T2)
    
    func=lambda x: (1/kappa(x))
    func_C=lambda x: (1/kappa_C(x))
    
    T_1.append(Factor*quad(func,T1,T2)[0])
    T_C1.append(Factor*quad(func_C,T1_C,T2_C)[0])

for i in range(len(f)):
    
    F=f[i]
    
    T1=Turning[i][1]
    T2=Turning[i][2]#
    
    T1_C=Turning_C[i][1]
    T2_C=Turning_C[i][2]#
    
    #print(T1,T2)
    
    func=lambda x: (1./kappa(x))
    func_C=lambda x: (1./kappa_C(x))
    
    T_2.append(Factor*quad(func,T1,T2)[0])
    T_C2.append(Factor*quad(func_C,T2,T2_C)[0])
#    T.append(Factor*quad(func,T1,T2)[0])
#    T_C.append(Factor*quad(func_C,T1_C,T2_C)[0])

plt.subplot(211)
plt.plot(f,T_1,"r",label="Region 1: uncorrected")
plt.plot(f,T_C1,"r--",label="Region 1: corrected")
plt.ylabel("Time [as]")
plt.title("Traversal time")
plt.legend()
plt.subplot(212)
plt.plot(f,T_2,"r",label="Region 2: uncorrected")
plt.plot(f,T_C2,"r--",label="Region 2: corrected")
plt.ylabel("Time [as]")
plt.xlabel("Field (a.u)")
plt.legend()
plt.savefig("traversal_time.png")
plt.show()
plt.close()



#
#W_alex=[8.73,10.19,11.96,13.55,15.76,17.22,19.08,20.85]
#T_alex=[34.30,39.83,46.29,51.85,60.15,65.68,73.04,79.50]
#
#plt.plot(W,f,label="Width with uncorrected")
#plt.plot(W_C,f,label="Width with corrected")
#plt.xlabel("Barrier width")
#plt.ylabel("Field")
#plt.legend()
#plt.show()
#plt.close()
#
#
#
#plt.plot(W,Time,label="Unconrrected")
#plt.plot(W_C,Time_C,label="Corrected")
#plt.plot(W_alex,T_alex,label="Alex")
#plt.xlabel("Barrier width")
#plt.ylabel("Time[as]")
#plt.legend()
#plt.show()
#plt.close()





















