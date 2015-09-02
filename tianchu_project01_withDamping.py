#Date: Jan 26 2015
#Author: Tianchu(Alex) Liang
#The following code implements double pendulum with damping.

#First, import necessary things:
import numpy as np
from scipy.integrate import odeint
from math import *

#Then, we define the following parameters for the system.
#The following parameters are measured by Daniel as actual values. 
m1=.339
m2=.339
l1=0.4325
l2=0.4325
g=9.81
k1=0.1
k2=0.1

#There are two versions of implementation
#that solve the equations of motion. The first one below
#simply defines df_dt normally. The second one after that 
#defines df_dt in dimensionless form. A derivation of dimensionless
#EOM is included in the manuscript. 
#====================
# Non-dimensionless Units
#====================
#==================================================
# A=(g*(m1/2+m2))/(l1*(m1/3+m2))
# B=(k1/(l1*(m1/3+m2)))
# C=(3*g)/(2*l2)
# D=3*k2/(m2*l2*l2)

# def df_dt(f,t):
# 	q1,q2,dq1,dq2=f
# 	ddq1=A*sin(q1)-B*dq1
# 	ddq2=C*sin(q2)-D*dq2
# 	return np.array([dq1,dq2,ddq1,ddq2])
#==================================================

#=====================
# Dimensionless:
# tau=sqrt((l1*(m1/3+m2))/(g*(m1/2+m2)))
#=====================
#==================================================
tau=sqrt((l1*(m1/3+m2))/(g*(m1/2+m2)))
A=k1/(l1*(m1/3+m2))
B=3*g/(2*l2)
C=3*k2/(m2*l2*l2)

def df_dt(f,t):
	q1,q2,dq1,dq2=f
	ddq1=sin(q1)-tau*A*dq1
	ddq2=B*tau*tau*sin(q2)-C*tau*dq2
	return np.array([dq1,dq2,ddq1,ddq2])
#==================================================

#Having defined df_dt, we define initial values for f and t.
#Notice, t here is dimensionless. The definition of dimensionless
#t corresponds to 20 secs in real time in 601 intervals. 
f0 = np.array([pi/18,0,0,0], dtype=np.float64)
t = np.linspace(0,20/tau,601, dtype=np.float64)

#To convert dimensionless t to actual time array, we multiply
#t by tau, and get tt.
tt = tau*t
#Solve the differential equations using df_dt, f0, tt:
f = odeint(df_dt, f0, tt)
#Write the result to a file:
outfile = open( "dbpend_damped.dat","w")
outfile.write( "Double Pendulum with damping\n")
outfile.write( "q1	q2	\n")
for q1, q2, dq1, dq2 in f:
  outfile.write( "{}	{}	\n".format(q1,q2))
outfile.close()