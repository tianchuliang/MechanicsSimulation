#Date: Feb7 2015
#Author: Tianchu(Alex) Liang
#IHRTLUHC

import numpy as np
from scipy.integrate import odeint

#Six Sets of Parameters:
a=0.1
b=0.1
c=14
# 
a1=0.1
b1=0.1
c1=4
# 
a2=a1
b2=b1
c2=6
# 
a3=a1
b3=b1
c3=8.5
# 
a4=a1
b4=b1
c4=8.7
# 
a5=a1
b5=b1
c5=9
# =======================
#Differential Equations:
def df_dt(f,t):
	x,y,z=f
	dx_dt=-y-z
	dy_dt=x+a5*y
	dz_dt=b5+z*(x-c5)
	return np.array([dx_dt, dy_dt, dz_dt])

#Define initial conditions of f and t,:
f0 = np.array([-10,0,0], dtype=np.float64)
t = np.linspace(0,1000,100000, dtype=np.float64)

#Solve the differential equations
f = odeint(df_dt, f0, t)

#Write the result to a file:
outfile = open( "R_attractor5.dat","w")
outfile.write( "Rossler Attractor\n")
outfile.write( "x y z\n")
for x,y,z in f:
  outfile.write( "{} {} {}\n".format(x, y, z))
outfile.close()