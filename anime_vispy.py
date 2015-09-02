import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from numpy import sin, cos, pi, array
from math import *


l1=1
l2=2

filename='/users/tianchu/documents/advmechcompphy/project1/dbpend.dat'
f=open(filename,'r')
clns=[]
for line in f:
    line=line.strip()
    column=line.split()
    clns.append(column)
clns.pop(0)
clns.pop(0)
x1=[]
x2=[]
y1=[]
y2=[]
for i in range(0,len(clns)):
	x1.append(sin(float(clns[i][0])))
	y1.append(cos(float(clns[i][0])))
	x2.append(sin(float(clns[i][0]))+sin(float(clns[i][1])))
	y2.append(cos(float(clns[i][0]))+cos(float(clns[i][1])))

x1=np.array(x1)
y1=np.array(y1)
x2=np.array(x2)
y2=np.array(y2)

print(len(x2))

fig = plt.figure()
ax = plt.axes(xlim=(-2.5, 2.5), ylim=(-2.5, 2.5))
ax.grid()
line1, = ax.plot([], [],'o-',lw=2)

def init():
    line1.set_data([],[])
    return line1,

def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]
    line1.set_data(thisx, thisy)
    return line1, 

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(x2), interval=1, blit=True)
plt.show()

