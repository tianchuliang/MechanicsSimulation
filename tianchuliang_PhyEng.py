#Author: Tianchu Liang
#Date: Feb22 2015
#IHRTLUHC
#
#-------------------------------------------------
#The program will iteratively find the roots where 
#collision happens. The ball drops from the initial position with
#initial velocity 0. Then, after the first collision, the ball change's its
#momentum. This leads to a new initial position and initial 
#velocity in both x and y direction. Then the program carries on looking
#for a new collision. 
#-------------------------------------------------
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import root,brentq
from scipy.integrate import odeint
interval=1/500
# Tfuture=10
delta=0.00001
xoriginal=0.0
yoriginal=10
voriginalx=0.
voriginaly=0.
g=9.8
# r=0.9 #Coefficient of Restitution
r=0.9 #Coefficient of Restitution
vx0=voriginalx #Initial velocity in x direction is 0; this variable maybe redefined later.
vy0=voriginaly #Initial velocity in y direction is 0; this variable maybe redefined later.
x0=xoriginal #Inditial x position; this variable maybe redefined later.
y0=yoriginal #Initial y position; this variable maybe redefined later.
#==================================================
#Define the surface function:
#Possible surface functs: 
#10-2*x+0.5*math.sin(8*x)
#10-4*x
#10-x*x
#(x-3)*(x-3)+1
#(1/5)*(x-5)*(x-5)
#(4/5)*(x-(5/2))*(x-(5/2))
def surface(x):
	#Make sure the surface function only goes to x=5.
	if(x>5):
		return 0
	else:
		#Make sure the surface function does not have negative altitude.
		if((4/5)*(x-(5/2))*(x-(5/2))<0):
			return 0
		else:
			return (4/5)*(x-(5/2))*(x-(5/2))
#==================================================
#Define the equation of motion for the bouncing ball:
def df_dt(f,t):
	x,y,vx,vy=f
	dx=vx
	dy=vy
	dvx=0 #Assuming no air resistance in the x direction
	dvy=-g #Assuming no air resistance in the y direction, only gravity.
	return (np.array([dx,dy,dvx,dvy]))
#==================================================
#Define the position of the bouncing ball at time t:
def info(t):
	f0 = np.array([x0,y0,vx0,vy0])
	times = np.array([0,t])
	f = odeint( df_dt, f0, times)
	x, y, vx, vy = f[-1]
	return( x, y,vx,vy)
#==================================================
#Define the y difference function: bouncing ball's y position minus the surface
#function at a certain location x. We look for a time t, when this funciton 
#is 0, namely a collision occurs.
# def y_diff(t_array):
# 	t = t_array[0]
# 	# print("time is {}".format(t))
# 	x, y, vx, vy= info(t)
# 	return( np.array([(surface(x))-y]))#The function surface is:1-x
def y_diff(t):
#	t = t_array[0]
	# print("time is {}".format(t))
	x, y, vx, vy= info(t)
#	return( np.array([(surface(x))-y]))#The function surface is:1-x
	# print("t={}, surface-y={}".format(t,y-(surface(x))))
	return(y-(surface(x)))#The function surface is:1-x


#==================================================
def find_tStep():
	if (vx0==0):
		vxThisMoment=0.001
	else:
		vxThisMoment=math.fabs(vx0)

	vyThisMoment=math.fabs(vy0)

	if(vyThisMoment>vxThisMoment):
		step=interval/vyThisMoment
	else:
		step=interval/vxThisMoment
	return (step)
#==================================================
#Define the impact function; it returns the x, y coordinate and time value for 
#the collision.
impactx=[]
impacty=[]
def impact():
	# sol = root(y_diff, [Tfuture])
	# t = float(sol.x)
	stop=0
	step=find_tStep()
	print("Step is: {}".format(step))
	tEnd=1.e-12
	while(stop==0):
		if(y_diff(tEnd)<0):
			print("!!!")
			stop=1
		else:
			tEnd=tEnd+step
	print("ydiff(1.e-12)={},ydiff(tEnd)={}".format(y_diff(1.e-12), y_diff(tEnd)))
	t = brentq(y_diff, 1.e-12, tEnd)
	x,y,vx,vy= info(t)
	impactx.append(x)
	impacty.append(y)
	return(x,y,vx,vy,t)
#==================================================
#Define the calculation for the angle between the slope 
#and the vertical line. Notice, a negative sign is necessary for calculation.
def cal_angle(x_pos):	
	yprev=surface(x_pos-delta)
	yafter=surface(x_pos+delta)
	rise=yafter-yprev
	run=2*delta
	return (-math.atan(rise/run))
#==================================================
#Transform current vx,vy into the correct next vx,vy, using vnormal,vtan:
def transform(vx,vy,angle):
	# print("angle is {}".format(angle))
	if(angle>0.):
		#vn is V in the normal direction: vx vy projected together to the normal direction of 
		#slope's tangent. 
		#vt is V in the tangent direction.
		vn=vx*math.sin(angle)+vy*math.cos(angle)
		vt=vx*math.cos(angle)-vy*math.sin(angle)
		#Adjust the value of vn,vt: 
		vt=vt
		vn=-r*vn
		#Transform vn,vt back to vx,vy:
		newvx=vn*math.sin(angle)+vt*math.cos(angle)
		newvy=vn*math.cos(angle)-vt*math.sin(angle)
		return(newvx,newvy)
	elif(angle<0.):
		vn=vx*math.sin(angle)+vy*math.cos(angle)
		vt=vx*math.cos(angle)-vy*math.sin(angle)
		#Adjust the value of vn,vt: 
		vt=vt
		vn=-r*vn
		#Transform vn,vt back to vx,vy:
		newvx=vn*math.sin(angle)+vt*math.cos(angle)
		newvy=vn*math.cos(angle)-vt*math.sin(angle)
		return(newvx,-newvy)
	else:
		vn=vx*math.sin(angle)+vy*math.cos(angle)
		vt=vx*math.cos(angle)-vy*math.sin(angle)
		#Adjust the value of vn,vt: 
		vt=vt
		vn=-r*vn
		#Transform vn,vt back to vx,vy:
		newvx=vn*math.sin(angle)+vt*math.cos(angle)
		newvy=vn*math.cos(angle)-vt*math.sin(angle)
		return(newvx,newvy)

#==================================================
# print(cal_angle(impact()[0]))
# print(impact()[2],impact()[3])
# print(transform(impact()[2],impact()[3],cal_angle(impact()[0])))
#==================================================
#This section defines the record() function. It takes in all the 
#necessary information to record the previous trajectory of the 
#bouncing ball as well as the surface. 
xpos=[]
ypos=[]
surfacePos=[]
def record(x,y,vx,vy,t):
	timeDuration=t
	fstart = np.array([x,y,vx,vy])
	timesteps=np.linspace(0,timeDuration,1000,dtype=np.float64)
	solu=odeint(df_dt,fstart,timesteps)
	for x,y,vx,vy in solu:
		xpos.append(x)
		ypos.append(y)
		surfacePos.append(surface(x))
	
#==================================================
#The following is kind of the "main" method in this Python program:
#It contains a while loop that finds/computes/records the first "counter"
#many collisions of the bouncing ball.
#I think looking for the first 20 collisions is sufficent.

counter=0#It counts how many collisions already happpend. We will restrict 
		 #our simulation to have at most 20 collisions.
while(counter<20):
	#Update the progress
	counter=counter+1
	col_x,col_y,col_vx,col_vy,col_t=impact()
	# print("{}th collision: x={},y={},vx={},vy={},t={}".format(counter,col_x,col_y,col_vx,col_vy,col_t))
	#Use the old initial condition, and the collision time, to record the trajectory
	#before the first collision:
	record(x0,y0,vx0,vy0,col_t)
	# print("Recorded trajectory starting from x={}, y={}, vx={}, vy={}, t=0 to {}".format(x0,y0,vx0,vy0,col_t))
	#Calculate the collision angle and Update the initial conditions:
	angle=cal_angle(col_x)
	# print("Collision angle is {}".format(angle))
	x0=col_x
	y0=col_y
	vx0,vy0=transform(col_vx,col_vy,angle)
	print("New info: {}, {}; {}, {}.".format(x0,y0,vx0,vy0))
	# print("New initial conditions: x0={},y0={},vx0={},vy0={}".format(x0,y0,vx0,vy0))


# Write out the data:
outputfile = open("pos.dat","w")
outputfile2 = open("col.dat","w")
outputfile3 = open("surface.dat","w")

#Put initial conditions into the first pos.dat file:
outputfile.write("Initial conditions of r x0 y0 vx0 vy0 are: \n")
outputfile.write("{} {} {} {} {} \n".format(r,xoriginal,yoriginal,voriginalx,voriginaly))

for x in range(0,len(xpos)):
  outputfile.write("{} {} \n ".format(xpos[x],ypos[x]))
for y in range(0,len(impactx)):
  outputfile2.write("{} {} \n".format(impactx[y],impacty[y]))
for z in range(0,len(xpos)):
  outputfile3.write("{} {} \n ".format(xpos[z],surfacePos[z]))


#Simple plotting section below for a quick peak:
plt.plot(xpos,ypos,'r',xpos,surfacePos,impactx,impacty,'bo')
plt.ylim([0,6])
plt.xlim([0,5])
plt.show()













