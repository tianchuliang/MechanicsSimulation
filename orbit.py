# Author: Jef Wagner
# Date: 20-02-2015

import numpy as np

# This file is an example file that shows you how to use scipy's root-
# finding function. It solves Kepler's equation, which is a
# trancendental equation and has no closed form answer. Therefore, we
# have to use a computer to solve the equation if we want a numerical
# value.  We will use the solution to Kepler's equation to find the
# positions along an orbit as a function of time.

# From the scipy package we load the `root` function.
from scipy.optimize import root

# Define the constants for the orbit
a = 2.5 # semi-major axis, 2.5 distance units
T = 10 # period, 10 time units (seconds I guess)
w = 2*np.pi/T # mean angular frequency of the orbit
e = 0.4 # eccentricity between, a number between 0 and 1

# The function that I'm trying to solve:
#
# E - e*sin(E) == w*t
#
# I'm trying to solve this equation for `E`, when `t` is given. I know
# (after reading scipy's documentation) that I need to manipulate this
# function so that it equals to zero when the equation is true. So I
# will define
#
# f(E,t) = E - e*sin(E) - w*t,
#
# Note that the variable I want to solve for (in this case `E`) has to 
# be the first variable. Furthermore the variable needs to be an array
# (this same `root` function generalizes to solving simultaneous 
# system of equations) and returns an array.
def f(E_array, t):
  E = E_array[0]
  return( np.array(E-e*np.sin(E)-w*t))

# In order to solve the above function I need to pass it to the
# `root` command. Below is an example of how that is done
t = 0.1 # the time will be 0.1
E0 = w*t # This is the initial guess
# We then call `root` scipy's root finding function.
# - The first argument is the function,
# - The second argument is the initial guess (as an array)
# - The third argument is other arguments of the `f` function
sol = root( f, [E0], t)
# The `root` function returns a "solution object", which I've 
# called `sol`. To get the value of our solution, I need to
# get the `x` member.
E = sol.x # Get the solution out of the solution object
print( "E(t={})={}".format(t,E)) # print it to the command line

# I want the previous steps to happen fairly automatically. This
# function will use the root command to find the eccentric anomaly,
# then use that to find the radius and angle (also called true
# anomaly).
def orbit(t):
  # We will use a constant angular speed as our initial guess
  E0 = w*t
  # Call the solution
  sol = root( f, [E0], t)
  # The value is in a numpy array, the float command just makes it
  # a normal number
  E = float(sol.x)
  # I then use the formula's to get `r` and `phi` and return them.
  r = a*(1-e*np.cos(E))
  phi = np.arctan2(np.sqrt(1-e**2)*np.sin(E),np.cos(E)-e)
  return( r, phi)

# Write out the orbital positions to a file:
outputfile = open( "orbit.dat", "w")
outputfile.write("Orbital Positions every 1/30th of a sec\n")
outputfile.write("t      r       phi\n")
# I use linspace to make 300 points between 0 and 10
for t in np.linspace(0,10,300,endpoint=False):
  r, phi = orbit(t)
  outputfile.write("{} {} {}\n".format(t,r,phi))
outputfile.close()