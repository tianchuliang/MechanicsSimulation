# Author: Jef Wagner
# Date: 20-02-2015

import numpy as np

# This file is an example file that shows how to string numerical
# methods together. In this example we try to find out at which 
# angles to fire a ping-pong ball in order to hit a target. Without
# air resistance this can be done analytically, but as soon as we
# include air restistance, we have to use a computer.

# We will use the general root finder, as well as the bracketing
# root finder  'brentq'.
from scipy.optimize import root, brentq
# We will also need the odeint to find the trajectories.
from scipy.integrate import odeint

# These are the parameters for our problem
a = 0.186 # 1/m (rho*cd*A)/(2*m) for ping pong ball
g = 10 # m/s^2 acceleration of gravity
v0 = 9 # m/s initial velocity
xf = 3 # m (x distance traveled)
yf = 1 # m (y target height)

# This is the differential equation for the ping pong ball
def dfdt( f, t):
  x,y,vx,vy = f
  dx = vx
  dy = vy
  # Drag is proportional to v^2 in the opposite direction of v
  dvx = -a*np.sqrt(vx*vx+vy*vy)*vx
  # and gravity is down
  dvy = -g -a*np.sqrt(vx*vx+vy*vy)*vy
  return( np.array([dx,dy,dvx,dvy]))

# This function finds the position of a ping-pong ball fired
# at `theta` after time `t`.
def pos( t, theta):
  vx0 = v0*np.cos(theta)
  vy0 = v0*np.sin(theta)
  f0 = np.array([0,0,vx0,vy0])
  times = np.array([0,t])
  f = odeint( dfdt, f0, times)
  x, y, vx, vy = f[-1]
  return( x, y)

# This is the first root finding function,
#
# f1(t, theta) = x(t)-xf
#
# It will return the `t` where f1(t,theta)=0.
def f1(t_array, theta):
  t = t_array[0]
  x, y = pos(t, theta)
  return( np.array([x-xf]))

# This function uses the root finding at what time x crosses xf, then
# finds the y value at that time.
def y_impact( theta):
  sol = root( f1, [1], theta)
  t = float(sol.x)
  x, y = pos( t, theta)
  return(y)

# Before solving for which theta give us a y_impact = yf, I decided
# to do a little exploration. This helped me figure out valid ranges
# in which to look for a solution. After I ran this, I wrote the
# bracketing values, then commented it out.
# for th in np.linspace(0,np.pi/2,5,endpoint=False):
#   y = y_impact( th)
#   print( "th={} y={}".format(th,y))

# This is the second root finding function,
#
# f2(theta) = y_impact(theta)-yf
#
# It will return the `theta` where f2(theta)=0.
def f2( theta):
  return( y_impact(theta)-yf)


# The backeting values are theta = 2pi/10 when ping_pong ball hits
# below the target and theta = 3pi/10 where it hits above.
th_lo = 2*np.pi/10
th_hi = 3*np.pi/10

# The  `brentq` command is for bracketed searches of 1-d functions. 
# - The first argument is the function to be minimized.
# - The second and third arguments are the bracketing values.
# The function simply returns the value of theta where f2=0.
th0 = brentq( f2, th_lo, th_hi)
# Print these values to the command line
y0 = y_impact(th0)
print( "th={}, y={}".format( th0, y0))

# Same thing but for now this time at 3 pi/10 the ping pong ball hits
# above and theta = 4pi/10 where it hits below.
th_lo = 3*np.pi/10
th_hi = 4*np.pi/10
th1 = brentq( f2, th_lo, th_hi)
y1 = y_impact(th1)
print( "th={}, y={}".format( th1, y1))
