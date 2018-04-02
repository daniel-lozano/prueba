import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from collections import Counter

name=argv[1]

F=open(name,"r")
text=F.readlines()
DATOS=[]
alfabeto=[]

for i in range(len(text)):
    for j in range(len(text[i])):
            DATOS.append(text[i][j])
    
print("datos a comprimir")
print(DATOS)

dictionary=Counter(DATOS)
print(dictionary.most_common(2))
print("len dic",len(dictionary))

array=dictionary.elements()
print(type(array))

"""
for i in range(1,len(array)):
    if(array[i]!=array[i-1]):
        alfabeto.append(array[i])
print(alfabeto)
"""
    

    
