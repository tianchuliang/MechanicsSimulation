#Date: Jan 24 2015
#Author: Tianchu(Alex) Liang

#First, import necessary things:
import numpy as np
from scipy.integrate import odeint
from math import *

#Then, we define the following parameters for the system.
#The following parameters are measured by Daniel as actual values. 
m=1
l=1
g=-9.8


#The one below
#simply defines df_dt:
def df_dt(f,t):
	q1,q2,dq1,dq2=f
	
	A=3*(9*g*sin(q1)+3*g*sin(q1-2*q2)+2*l*(2*dq2**2+3*dq1**2*cos(q1-q2))*sin(q1-q2))
	B=l*(-23+9*cos(2*q1-2*q2))
	C=-3*(16*l*dq1**2*sin(q1-q2)+3*l*dq2**2*sin(2*q1-2*q2)+9*g*sin(2*q1-q2)-7*g*sin(q2))

	ddq1=A/B
	ddq2=C/B

	return np.array([dq1,dq2,ddq1,ddq2])
#Having defined df_dt, we define initial values for f and t.
f0 = np.array([pi/18,0,0,0], dtype=np.float64)
t = np.linspace(0,50,5000, dtype=np.float64)
#Solve the differential equations using odeint()
f = odeint(df_dt, f0, t)
#Write data to a file:
outfile = open( "dbpend.dat","w")
outfile.write( "Double Pendulum\n")
outfile.write( "q1	q2	\n")
for q1, q2, dq1, dq2 in f:
  outfile.write( "{}	{}	\n".format(q1,q2))
outfile.close()