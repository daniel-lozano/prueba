import numpy as np
import matplotlib.pyplot as plt
from sys import argv

name_alfabetos=argv[1]
F=open(name_alfabetos,"r")
text=F.readlines()

alfabeto=[]
alfabeto_binario=[]

#-----------------------------Leyendo los alfabetos------------------------------
index=-1

for i in range(len(text)):
    
    
    if(text[i][0]!="\n" and i!=index):
        alfabeto.append(text[i][0])
        alfabeto_binario.append(text[i][2:-1])
    
    if(text[i][0]=="\n" and index==-1):
        alfabeto.append(text[i][0])
        alfabeto_binario.append(text[i+1][1:-1])
        index=i+1
#--------------------------Leyendo texto codificado----------------------------------
name_code=argv[2]
F=open(name_code,"r")
text=F.readlines()
DATOS=[]

for i in range(len(text)):
    for j in range(len(text[i])):
        DATOS.append(text[i][j])
F.close()

#-------------------------Decodificnaod------------------------------------------
name="decode.txt"
F=open(name,"w")
palabra=""
for i in range(len(DATOS)):
    
    if(DATOS[i]=="0"):
        palabra+=DATOS[i]
    if(DATOS[i]=="1"):
        palabra+=DATOS[i]
        pos=alfabeto_binario.index(palabra)
        F.write(alfabeto[pos])
        palabra=""
F.close()



















