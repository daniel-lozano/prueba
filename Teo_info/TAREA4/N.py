import numpy as np
import matplotlib.pyplot as plt
from sys import argv

N=int(argv[1])

def function(N):
    f=np.zeros(2**N+1)
    dim=np.zeros(2**N+1)
    for i in range(len(f)):
        f[i]=np.log(2**N-i)/(N*np.log(2.0))
        dim[i]=i/2.0**N
    return f,dim

for i in range(5):
    
    F,D=function(N+i)
    plt.plot(F,D,".-",label=str(N+i))

plt.legend()    
plt.show()
