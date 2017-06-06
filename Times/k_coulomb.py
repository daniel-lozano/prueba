import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

hbar=197#MeV
m=0.51 #MeV

w=1# some value
F=1# some value

t=1 #some value
E=0.5#some value



def func(x,t,E,m,hbar,F,w):
    x0=0.025
    
    fact1=pow(2*m,0.5)/hbar
    
    func1=-1/np.sqrt(x**2+x0**2)
    
    func2=-x*F*np.cos(w*t)-E
    
    return fact1*np.sqrt(abs(func1+func2))

def turning(E,m,F,t,w):

    x1=1
    x2=2
    arreglo=np.zeros(2)
    arreglo=[x1,x2]
    return arreglo


#------------corriendo el programa----------------#



X=turning(E,m,F,t,w)
funcion= lambda x: -2*func(x,t,E,m,hbar,F,w)


resultados=integrate.quad(funcion,X[0],X[1],limit=100)
