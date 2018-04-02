import numpy as np
import matplotlib.pyplot as plt
from sys import argv

name=argv[1]

F=open(name,"r")
text=F.readlines()
DATOS=[]

for i in range(len(text)):
    for j in range(len(text[i])):
            DATOS.append(text[i][j])
            
    
print(DATOS)


    

    
