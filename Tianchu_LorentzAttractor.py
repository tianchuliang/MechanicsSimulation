#Date: Jan23 2015
#Author: Tianchu(Alex) Liang
#IHRTLUHC

import numpy as np
from scipy.integrate import odeint

#Parameters:
sigma=50
beta=5
rho=80

#Differential Equations:
def df_dt(f,t):
	x,y,z=f
	dx_dt=(sigma*(y-x))
	dy_dt=(x*(rho-z)-y)
	dz_dt=(x*y-beta*z)
	return np.array([dx_dt, dy_dt, dz_dt])

#Define initial conditions of f and t:
f0 = np.array([0.1,1.2,0.5], dtype=np.float64)
t = np.linspace(0,8.15,100000, dtype=np.float64)

#Solve the differential equations
f = odeint(df_dt, f0, t)

#Write the result to a file:
outfile = open( "attractor.dat","w")
outfile.write( "Lorentz Attractor\n")
outfile.write( "x y z\n")
for x,y,z in f:
  outfile.write( "{} {} {}\n".format(x, y, z))
outfile.close()





