import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from collections import Counter

name=argv[1]
F=open(name,"r")
text=F.readlines()
DATOS=[]
#-----------------------------Leyendo el texto------------------------------

for i in range(len(text)):
    for j in range(len(text[i])):
            DATOS.append(text[i][j])
F.close()
#-----------------------------Generando un diccionario y alfabetos de codigo plano------------------------------
alfabeto=[]
frec_alf=[]
binario=[]

dictionary=Counter(DATOS)
array=dictionary.most_common()

def copying_but_last(A):
    B=[]
    for i in range(len(A)-1):
        B.append(A[i])
    return B

def searching(A,a):
    
    for i in range(len(A)):
        if(A[i]==a):
            return i
    return -1




letras=[]
valores=[]
for i in range(len(array)):
    letras.append(array[i][0])
    valores.append(array[i][1])
    
print(letras,valores)







"""
array2=copying_but_last(array)
print(array)
print(array2)
"""


