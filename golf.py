#Author: Tianchu Liang	
#Date: 1/3/2015

import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d,UnivariateSpline
import matplotlib.pyplot as plt
from scipy.optimize import root, brentq, minimize
import math
#================Choose problem:==============================
Golf_Or_Ping=0# 0 for golfball, 1 for ping pong ball.
dimpled_Or_Smooth=0# 0 for dimpled, 1 for smooth.
#================Define Constants:============================
m_club=0.2#kg
v0_club=40#Initial velocity of the golf club, in m/s
rho=1.275#kg/m^3
mu=1.845e-5 #kg*sec^-1*m^-1
g=9.8
t_height=0.02#meters
groundLevel=0#meter
#------------------------
m_golf=0.04593#kg
m_ping=0.0027#kg
#------------------------
D_golf=0.04267#meter
D_ping=0.04#meter
#------------------------
A_golf=math.pi*(D_golf/2)**2
A_ping=math.pi*(D_ping/2)**2
#------------------------
#Assuming elastic collision between the golf club and the golf ball, 
#we can derive the initial velocity of golf ball to be v0=(m_club/m_golf)*v0_club
# v0=
v0=2*m_club*v0_club/(m_club+m_golf)
#------------------------
dimplex=[]
dimpley=[]
smoothx=[]
smoothy=[]
xtraject=[]
ytraject=[]
#====================Golf Ball or Ping Pong Ball===============
m=0
D=0
A=0
if(Golf_Or_Ping==0):
  #Golf:
  print("Golf")
  m=m_golf
  D=D_golf
  A=A_golf
else:
  print("PingPong")
  #PingPong:
  m=m_ping
  D=D_ping
  A=A_ping
#====================Golf Ball or Ping Pong Ball===============
#====================Load up the data:=========================
filename1='/users/tianchu/documents/advmechcompphy/project2/dimpled.dat'
filename2='/users/tianchu/documents/advmechcompphy/project2/smooth.dat'
f1=open(filename1,'r')
f2=open(filename2,'r')
clns1=[]
clns2=[]

for line in f1:
    line=line.strip()
    column=line.split()
    clns1.append(column)

for line in f2:
    line=line.strip()
    column=line.split()
    clns2.append(column)

clns1.pop(0)
clns2.pop(0)

for i in range(0,len(clns1)):
	dimplex.append(float(clns1[i][0]))
	dimpley.append(float(clns1[i][1]))

for j in range(0,len(clns2)):
	smoothx.append(float(clns2[j][0]))
	smoothy.append(float(clns2[j][1]))

#Modify dimplex, dimpley and smoothx, smoothy in order to
#make UnivariateSpline work well.
insert1=dimplex[0]
insert2=dimpley[0]
insert3=smoothx[0]
insert4=smoothy[0]
append1=dimplex[-1]
append2=dimpley[-1]
append3=smoothx[-1]
append4=smoothy[-1]
dimplex.insert(0,insert1)
dimpley.insert(0,insert2)
smoothx.insert(0,insert3)
smoothy.insert(0,insert4)
dimplex.insert(-1,append1)
dimpley.insert(-1,append2)
smoothx.insert(-1,append3)
smoothy.insert(-1,append4)
#===============Interpolate functions:=========================
f_dimple=UnivariateSpline(dimplex,dimpley)
f_smooth=UnivariateSpline(smoothx,smoothy)
#-------Choose the interpolated function:---------
if(dimpled_Or_Smooth==0):
  f_inter=f_dimple
else:
  f_inter=f_smooth
#-------------------------------------------------
#==============================================================
def df_dt( f, t):
  x,y,vx,vy = f

  Re=(rho*np.sqrt(vx**2+vy**2)*D)/mu

  dx = vx
  dy = vy

  # print(Re)
  
  c_d=f_inter(Re)

  a_dragx=(1/2)*c_d*(rho/m)*A*np.sqrt(vx*vx+vy*vy)*vx
  a_dragy=(1/2)*c_d*(rho/m)*A*np.sqrt(vx*vx+vy*vy)*vy

  # Drag is proportional to v^2 in the opposite direction of v
  dvx = -a_dragx
  # and gravity is down
  dvy = -g - a_dragy

  return(np.array([dx,dy,dvx,dvy]))
#==============================================================
def pos(t,theta):
  vx0 = v0*np.cos(theta)
  vy0 = v0*np.sin(theta)
  f0 = np.array([0,t_height,vx0,vy0])
  times = np.array([0,t])
  f = odeint( df_dt, f0, times)
  x, y, vx, vy = f[-1]
  return( x, y)
#==============================================================
def hit_ground(t_array, theta):
  t = t_array[0]
  x, y = pos(t, theta)
  return( np.array([y-groundLevel]))
#==============================================================
def how_far(theta):
  sol = root( hit_ground, [10.0], theta)
  t = float(sol.x)
  x, y = pos( t, theta)
  return(x)
#==============================================================

for th in np.linspace(math.pi/10,2*math.pi/5,10,endpoint=False):
  x = how_far(th)
  print( "theta={} x={}".format(th,x)) 



