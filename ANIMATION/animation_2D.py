"""
    =================
    An animated image
    =================
    
    This example demonstrates how to animate an image.
    """
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()


def f(x, y, t):
    return (np.sin(x)*np.cos(y))*np.exp(-t**2/10)

x = np.linspace(0, 2 * np.pi, 100)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
t=0
vx=0.5
vy=-0.1

im = plt.imshow(f(x, y,t), animated=True)


def updatefig(*args):
    global x, y, t
    t+=1./50
    x+=vx*t
    y+=vy*t
    im.set_array(f(x, y,t))
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=50)
plt.xlabel("$ x $")
plt.ylabel("$ y $")

plt.show()

