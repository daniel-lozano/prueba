import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt





def inte_exp(x):
    res=[]
    func=lambda y: 1
    for i in range(len(x)):
        res.append(np.exp(-2*quad(func,0,x[i])[0]))
    return res


def inte_num(x):
   
    func=lambda y: 1
    return (np.exp(quad(func,0,x)[0]))
   

x=np.linspace(0,1)
y=inte_exp(x)
plt.plot(x,y,label="exp")
plt.show()
plt.close()
funcion=lambda z: inte_num(z)
expo=lambda z: np.exp(z)

t=np.linspace(0,1)
ft=[]
ft_prueba=[]


print("quad result=",quad(funcion,0,1)[0])
print("quad result=",quad(expo,0,1)[0])

for i in range(len(t)):
    ft.append(quad(funcion,0,t[i])[0])
    ft_prueba.append(np.exp(t[i])-1)
    

plt.plot(t,ft,label="function")
plt.plot(t,ft_prueba,label="desired result")
plt.legend()
plt.show()
plt.close()

    
