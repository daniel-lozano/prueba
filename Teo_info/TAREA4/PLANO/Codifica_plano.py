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
array=dictionary.most_common() #arreglo que contiene el diccionario del texto
numbers=np.zeros(len(array))

CODIGO=open("codigo_plano_"+name,"w")

for i in range(len(array)):
    #Encontrando el alfabeto del texto
    alfabeto.append(array[i][0])
    
    #Hallando las frecuencias de cada caracter
    frec_alf.append(array[i][1]*1.0/len(array))
    
    #Asignando valores para el histograma
    numbers[i]=i
    
    #Creando un alfabeto binario
    caracter=""
    for j in range(i):
        caracter+="0"
    caracter+="1"
    binario.append(caracter)
    CODIGO.write(alfabeto[-1]+" "+caracter+"\n")

CODIGO.close()


#---------------------------Generando histograma--------------------------
'''
datos_hist=np.ones(len(DATOS))

for i in range(len(DATOS)):
    
    for j in range(len(alfabeto)):
        
        if(DATOS[i]==alfabeto[j]):
            
            datos_hist[i]=numbers[j]

plt.hist(datos_hist,bins=len(array))
plt.xticks(numbers,alfabeto)
plt.show()
'''

#----------------------Escribiendo texto codificador---------------------------

FILE=open("comp_plano_"+name,"w")
for i in range(len(DATOS)):
    index=alfabeto.index(DATOS[i])
    palabra=binario[index]
    FILE.write(palabra)
FILE.close()


#---------------------Decodificando el texto------------------------------------




