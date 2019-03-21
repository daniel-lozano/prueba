
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import animatplot as amp

def D_func(X,Y,T):
    vx=0.1
    vy=-0.2
    return (np.cos(2*(X+vx*T))*np.sin(3*(Y+vy*T)))*np.exp((X+vx*T)**2/10)*np.exp(-T/10)


x = np.linspace(0, 1, 50)
t = np.linspace(0, 1, 50)

X, T = np.meshgrid(x, t)
Y = np.sin(2*np.pi*(X+T))

#Passing the image as a block to amp
block = amp.blocks.Line(X, Y)
#Animating the block
anim = amp.Animation([block])

#Add a pause botton and a bar to manage the time
anim.controls()
#anim.save_gif('line_moving') # save animation for docs
plt.show()

#Making a 2D animated plot with an evoolution variable t

x = np.linspace(-2, 2, 41)
y = np.linspace(-2, 2, 41)
t = np.linspace(0, 2*np.pi, 50)

X, Y, T = np.meshgrid(x, y, t)

pcolormesh_data =D_func(X,Y,T)
pcolormesh_block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data,
                                         t_axis=2, vmin=-1, vmax=1)
plt.colorbar(pcolormesh_block.quad)
timeline = amp.Timeline(t, fps=10)

# now to contruct the animation
anim = amp.Animation([pcolormesh_block], timeline)
anim.controls()
plt.show()






