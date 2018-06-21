import numpy as np
import matplotlib.pyplot as plt

FILE=open("data.dat","r")
Lines=FILE.readlines()
x=[]
y=[]
for line in Lines:
    x.append(float(line.split(",")[0]))
    y.append(float(line.split(",")[1].split("\n")[0]))
    

plt.scatter(x,y,label="Factorial")
plt.legend()
plt.savefig("Factorial.png")
