from sys import argv


a=int(argv[-1])
FILE=open("data.dat","a+")
fac=1
for i in range(1,a+1):
    fac*=i

print(argv[-1]+"!="+str(fac))
FILE.write(argv[-1]+","+str(fac)+"\n")
