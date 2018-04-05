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
print array
