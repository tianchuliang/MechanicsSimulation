# Author: Jef Wagner
# Date: 20-02-2015

import numpy as np

# This file is an example file that shows how to use scipy's
# optimization function. It finds the equilbrium configuration for a
# chain of masses with springs strung between them. After a long
# enough time, the system will find a settle into a minimum energy
# state, so we will find the configuration where the enery is a
# minimum.

# From the scipy package we loat the `minimize` function.
from scipy.optimize import minimize

# Define the constants for the system
l0 = 1. # unstretched length of the springs
m = 0.1 # value of all masses
g = 10. # gravity
k = 10. # spring constant

# I need to know a little about the function the structure of a
# function that I am going to minimize. The function has the form
#
# f( array, option_arg0, optional_arg1, ...):
#
# The first argument needs to be a 1-d array. It is with respect to
# this argument that the function will be minimized.

# The function that I want to minimize is the potential energy of the
# chain. It is simply the sum of the gravitational potential energies
# and the spring potential energies. These depend only on the (x,y)
# position of each mass, and the position of the two endpoints of the
# chain (which we will treat as fixed). Such a function is easy enough
# to write.
def chain_energy( pos, initial_pos, final_pos):
  # the total energy is a sum, so we start the sum at 0
  U = 0.
  # We add the gravitational potential energy for each mass.
  for x, y in pos:
    U += m*g*y
  # We add the first spring's potential energy
  dx, dy = pos[0]-initial_pos
  l = np.sqrt(dx*dx+dy*dy)
  U += 1/2*k*(l-l0)**2
  # We add the potential energy for each spring between the massess
  for i in np.arange(1,len(pos)):
    dx, dy = pos[i]-pos[i-1]
    l = np.sqrt(dx*dx+dy*dy)
    U += 1/2*k*(l-l0)**2
  # We add the potential energy of the last spring.
  dx, dy = final_pos-pos[-1]
  l = np.sqrt(dx*dx+dy*dy)
  U += 1/2*k*(l-l0)**2
  return( U)

# However, this function is not quite in the form that is needed by
# scipy's `minimize` function: the first argument in the above
# function is a 2-d array, but it should be a 1-d array. I can get
# arround that writing a very small wrapper function that changes
# the shape of the array using numpy's reshape method. The following
# function takes a 1-d array, then reshapes to a 2-d array.
def min_func( pos_array, *args):
  size = len(pos_array)/2
  return( chain_energy( pos_array.reshape(size,2), *args))


# In order to numerically minimize the system it must be well defined.
n = 20 # I will minimize the shape with 20 springs (19 masses)
xl = 0.; yl = 0. # The left hand endpoint will be at the origin
pos_l = np.array([xl,yl])
xr = 10.; yr = 0. # The right hand endpoint will be at (10,0)
pos_r = np.array([xr,yr])

# For my initial guess, I will simply guess that the masses will
# be spaced evenly between the endpoints.
x_pos = np.linspace(xl, xr, n, endpoint=False)[1:]
y_pos = np.linspace(yl, yr, n, endpoint=False)[1:]
# Note: the `[1:]` after the linspace command is called a "slice". You
# can find a lot of information about python slices online. The
# command as I used it drops the first element of the sequence, (which
# in this case would be the point x=0, y=0).
pos0 = np.array([[x,y] for x,y in  zip(x_pos, y_pos)])
# This command weaves the x,y positions together using something
# called a "list comprehension". Again, you can find plenty of
# information about python list comprehension online.

# We then call the `minimize` function:
# - The first argument is the function to be minimized
# - The second argument is the initial guess. Note: it should be a 1-d
#   array, so I use numpy's flatten method.
# - The third argument is a tuple of additional arguments for the 
#   function to be minimized
sol = minimize( min_func, pos0.flatten(), (pos_l, pos_r))
# The `minimize` function returns a "solution object", which I've
# called `sol`. To get the actual values I need to get the `x` member.
pos = sol.x
# Now `pos` has a 1-d array of x and y values. To make a 2-d array
# of (x,y) pairs, I can use numpy's reshape method.
size = len(pos)/2
pos = pos.reshape(size,2)

# Write out the positions of the endpoints and masses.
outputfile = open("spring_chain.dat","w")
outputfile.write("Position of first endpoint:\n")
outputfile.write("{} {}\n".format(pos_l[0], pos_l[1]))
outputfile.write("Number of masses: {}\n".format( n-1))
outputfile.write("Position of masses:\n")
for x, y in pos:
  outputfile.write("{} {}\n".format(x,y))
outputfile.write("Position of last endpoint:\n")
outputfile.write("{} {}\n".format(pos_r[0], pos_r[1]))
outputfile.close()