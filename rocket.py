# Author: Jef Wagner
# Date: 22-01-2015

# This file demonstrates the use of a ordinary differential equation
# solver. We will use it to plot the flight of a rocket in x-y space,
# where I treat the variable trust curve, the thrust acting along the
# body of the rocket, drag due to air resistance, and gravity. This
# involves solving Newtons equations of motion, which can be written
# as a system of coupled non-linear second order ordinary differential
# equations (ODE), with initial conditions.

# We have to import the NumPy and SciPy packages.
import numpy as np
from scipy.integrate import odeint

# I will pretend that my rocket has a frontal area of about 4 square
# centimeters
A = 4.0e-4
# A coefficent of drag of 0.3 (standard for small rockets)
c_d = 0.3
# A mass of 100 grams
m = 0.1
# And gravity is 10 meters/second squared
g = 10.0

# I modeled this thrust curve off of an Estes B4 model rocket motor.
# The thrust curve has a peak of around 12 Newtons shortly after
# takeoff (about 0.1 seconds) then is steady at about 3 newtons for 1
# second. The thrust is aligned along the body of the rocket, which is
# given by the direction of motion.
def F_thrust(t, vx, vy):
  if t < 0.2:
    F = 1200*t*(0.2-t)+3.
  elif t < 1.0:
    F = 3.
  else:
    F = 0.0
  v_mag = np.sqrt(vx**2+vy**2)
  return (F*vx/v_mag, F*vy/v_mag)

# I will treat the drag as simply aerodynamic drag, that is
# proportional to the square of the velocity.
def F_drag(vx, vy):
  v_mag = np.sqrt(vx**2+vy**2)
  pf = -0.5*c_d*A*v_mag
  return (pf*vx, pf*vy)

# The total force is simply the sum of the thrust and drag forces and
# gravity which only acts in the negative y direction.
def F_tot(t, vx, vy):
  Ftx, Fty = F_thrust(t, vx, vy)
  Fdx, Fdy = F_drag(vx, vy)
  Fx = Ftx+Fdx
  Fy = Fty+Fdy-m*g
  return (Fx, Fy)

# Now we have set up our problem, so we need to start thinking about
# what we need for the ODE solver `odeint`. We know that we have to
# feed the solver a callable function, and the initial condition, and
# we can see from the documentation, that we also have to feed in a
# list or array of times (the independent variable). The output fill
# be an array of values the same length as the list of times,
# corresponding to each time.

# First we define the callable function. The function only accepts a
# system of first order differential equations. The bulk of the
# conceptual difficulty comes here, I have to turn the second order
# differential equation into a system of first order differential
# equations. This is done by treating the positions `x` and `y` as
# different from the velocities `vx` and `vy`. We then need to define
# a function that takes the vector (vx, vy, x, y) and t as arguments,
# and returns a vector of the derivatives. The following function does
# exactly that for us.
def df_dt(f,t):
  vx, vy, x, y = f
  Fx, Fy = F_tot( t, vx, vy)
  dvx_dt = Fx/m
  dvy_dt = Fy/m
  dx_dt = vx
  dy_dt = vy
  return np.array([dvx_dt, dvy_dt, dx_dt, dy_dt])

# This is the array of initial values, I cannot set the initial
# velocity to be zero for two reasons, one practical, the other
# physical. Physically, I need the rocket to be moving in order to the
# velocity to have a definite direction. Practically, I divide by the
# magnitude of the velocity, so I will run into problems if the
# velocity is zero. So I set the initial velocity vector to be small,
# and slightly tilted off vertical.
f0 = np.array([0.01*np.sin(0.05),0.01*np.cos(0.05),0,0], dtype=np.float64)

# We will write 200 points between 0 and 8.15 seconds. I choose 8.15
# seconds, by running the program initially with 10, but the rocket
# dropped below y=0, and I only wanted to run until it hit the ground
# again.
t = np.linspace(0,8.15,200, dtype=np.float64)

# This is the command the preforms the numerical integration.
f = odeint( df_dt, f0, t)

# Then I write the x and y data to a file so that I can plot it.
outfile = open( "rocket.dat","w")
outfile.write( "Rocket Trajectory\n")
outfile.write( "x       y\n")
for vx, vy, x, y in f:
  outfile.write( "{} {}\n".format(x, y))
outfile.close()
