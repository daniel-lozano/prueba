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

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

R=1
g=9
w=2
m=1
n=(R*w**2)/g

delta = 0.025
#theta
x = np.arange(-1.1*np.pi, 1.1*np.pi, delta)
#p_theta
y = np.arange(-10.0, 10.0, delta)

X, Y = np.meshgrid(x, y)


#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
Z = (Y**2)/(2*m*R**2)+m*g*R*( 1-np.cos(X) - (n*np.sin(x)**2)/2  )


# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
fig=plt.figure()

levels = np.arange(-5.0,100,1)
CS = plt.contour(X, Y, Z, levels=levels)
#plt.clabel(CS, inline=1, fontsize=10)
plt.title("$ \eta < 1 $")

plt.show()
fig.savefig("eta<1.png")
#plt.savefig("eta<1.png")


w=7
n=(R*w**2)/g

delta = 0.025
#theta
x = np.arange(-1.1*np.pi, 1.1*np.pi, delta)
#p_theta
y = np.arange(-10.0, 10.0, delta)

X, Y = np.meshgrid(x, y)


#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
Z = (Y**2)/(2*m*R**2)+m*g*R*( 1-np.cos(X) - (n*np.sin(x)**2)/2  )


# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
fig=plt.figure()





CS = plt.contour(X, Y, Z, levels=levels)

#plt.clabel(CS, inline=1, fontsize=10)
plt.title("$ \eta > 1 $")

plt.show()
fig.savefig("eta>1.png")


