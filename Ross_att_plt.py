import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

mpl.rcParams['legend.fontsize'] = 10


filename='/users/tianchu/documents/advmechcompphy/MidTermTakeHome/Data/R_attractor5.txt'
f=open(filename,'r')
clns=[]
for line in f:
    line=line.strip()
    column=line.split()
    clns.append(column)

clns.pop(0)
clns.pop(0)


x=[]
y=[]
z=[]
for i in range(0,len(clns)):
	x.append(float(clns[i][0]))
	y.append(float(clns[i][1]))
	z.append(float(clns[i][2]))

# x=np.array(x)
# y=np.array(y)
# z=np.array(z)

print(type(x[1]))
print(len(x))


fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(x, y, z, c=u'y',label='Rossler Attractor')
ax.legend()
plt.show()