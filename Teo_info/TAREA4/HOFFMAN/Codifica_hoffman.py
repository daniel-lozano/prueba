import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from collections import Counter
from copy import deepcopy

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

letras=[]
valores=[]
codigo=[]

for i in range(len(array)):
    letras.append(array[i][0])
    valores.append(array[i][1])
    codigo.append("")

#-------------------------Useful functions------------------------------------------------------
def searching(A,a):
    for i in range(len(A)):
        if(A[i]==a):
            return True
    return False

def min_2(A):
    min1=min(A)
    index1=np.argmin(A)
    A[index1]=np.inf
    min2=min(A)
    index2=np.argmin(A)
    return [index1,index2]

def reorganiza(A):
    C=[]
    marker=False
    for i in range(len(A)-1):
        if(A[i][1]>=A[-1][1] or marker==True):
            C.append(A[i])
        
        if(A[i][1]<A[-1][1] and marker==False):
            C.append(A[-1])
            C.append(A[i])
            marker=True
    return C

def repite(A):
    
    for i in A:
        cont=0
        for j in A:
            if(i==j):
                cont+=1
        if(cont>=2):
            return True
    return False
#---------------------------preparing the new alphabet-------------------------------------------
old=[]
for i in range(len(array)):
    old.append([array[i][0],array[i][1]])

#---------------------------Writing the Binary code--------------------------------------------

for i in range(len(array)-1):
    
    for j in range(len(letras)):
        
        if(searching(old[-2][0],letras[j])==True):
            codigo[j]+="1"
        if(searching(old[-1][0],letras[j])==True):
            codigo[j]+="0"
    
    
    act=old[:-1]
    act[-1][0]+=old[-1][0]
    act[-1][1]+=old[-1][1]
    
    old=reorganiza(act)
    new=deepcopy(old)
print(letras)
print(codigo)
print(repite(codigo))

#-----------------------Writing code alphabet-------------------------------------------------

CODIGO=open("codigo_"+name,"w")
for i in range(len(letras)):
    CODIGO.write(letras[i]+" "+codigo[i]+"\n")
CODIGO.close()

COMP=open("comp_"+name,"w")

for i in range(len(DATOS)):
    for j in range(len(letras)):
        if(DATOS[i]==letras[j]):
            COMP.write(codigo[j])

COMP.close()




