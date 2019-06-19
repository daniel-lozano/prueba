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
    
    ro=1.0/(2*Z)
    
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
num=9

F=f[num]
print(F)
T1=Turning[num][1]
T2=Turning[num][2]
T1_C=Turning_C[num][1]
T2_C=Turning_C[num][2]

n=np.linspace(0.1,T2+2,1000)

V_ori=potential(n)
V_ori_triangle=potential_triangle(n,T1,T2)
V_cor=potential_schro(n)
V_cor_triangle=potential_triangle_C(n,T1_C,T2_C)

index1=np.where(n<T1)
index2=np.where(n<T1_C)
V_ori_triangle[index1]=0
V_cor_triangle[index2]=0


plt.plot(n[0:1000:10],V_cor[0:1000:10],".-",label="Corrected")
plt.plot(n,V_ori,":",label="Original",linewidth=2)

plt.plot(n,V_ori_triangle,"--",label="Original triangle")
plt.plot(n,V_cor_triangle,"--",label="Corrected triangle")
plt.legend()
plt.ylim(-0.3,0.3)
plt.xlim(min(n),15)
plt.xlabel("$ \eta $", size=15)
plt.ylabel("$ V(\eta,F) $", size=15)
plt.title("$F=0.11 $ ")
plt.savefig("triangle_potential_"+str(F)+".eps")

plt.show()



















