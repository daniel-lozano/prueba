"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

R=1
g=9
w=2
m=1
n=(R*w**2)/g

delta = 0.005
#theta
x = np.arange(-2 , 2 , delta)
#p_theta
y = np.arange(-2 , 2 , delta)

X, Y = np.meshgrid(x, y)

for i in range(4):

    Z =-1/np.sqrt(X**2+Y**2+0.025)-X*np.cos(i*np.pi/4)

    fig=plt.figure(1)

    levels = np.arange(-10,2,0.1)
    CS = plt.contour(X, Y, Z, levels=levels)

    plt.xlim(-2 ,2 )
    plt.ylim(-2 ,2 )
    plt.xlabel("$ x $")
    plt.ylabel("$ y $")
    plt.title("$ Potential\ contour\ plot\ $")

    if(i==0):
        plt.savefig("contour_potential.png")
    
    plt.show()
    plt.close()
   

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_wireframe(X, Y, Z, rstride=30, cstride=30)
    ax.set_xlabel("$ x $")
    ax.set_ylabel("$ y $")
    ax.set_zlabel("$ V(x,y,t) $")
   

    if(i==0):
       plt.savefig("3D_potential.png")
    plt.show()
