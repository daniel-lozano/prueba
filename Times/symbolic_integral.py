from sympy import *
import numpy as np

r=Symbol("r")
r0=Symbol("r0")
e=Symbol("e")
Ip=Symbol("Ip")
r=Symbol("r")


print(integrate(exp(r)))

result=integrate(2*r0*(r**4-(r0*r)**2 + ((e**2)*r**3/Ip)*(1+(1+r/(2*r0))*exp(-r/r0)))**(-0.5),r)

print( result)#,(r,r0,oo))
