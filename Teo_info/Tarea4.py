import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from collections import Counter

name=argv[1]

F=open(name,"r")
text=F.readlines()
DATOS=[]
alfabeto=[]
frec_alf=[]
binario=[]
#-----------------------------Leyendo el texto------------------------------

for i in range(len(text)):
    for j in range(len(text[i])):
            DATOS.append(text[i][j])
    
#-----------------------------Generando un diccionario y alfabeto------------------------------

dictionary=Counter(DATOS)

print("longitud del alfabeto",len(dictionary.most_common()))


array=dictionary.most_common()
numbers=np.zeros(len(array))
for i in range(len(array)):
    alfabeto.append(array[i][0])
    frec_alf.append(array[i][1]*1.0/len(array))
    numbers[i]=i

#---------------------------Generando histograma--------------------------
datos_hist=np.ones(len(DATOS))

for i in range(len(DATOS)):
    
    for j in range(len(alfabeto)):
        
        if(DATOS[i]==alfabeto[j]):
            
            datos_hist[i]=numbers[j]

plt.hist(datos_hist,bins=len(array))
plt.xticks(numbers,alfabeto)



