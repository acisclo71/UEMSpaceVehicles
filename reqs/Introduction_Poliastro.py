#!/usr/bin/env python
# coding: utf-8

# # Introduction to Python in Astrodynamics
# ## Poliastro

# For this first example, we are going to use the Orbit objects inside the poliastro.twobody module. They store all the required information to define and orit.
# * The body acting as the central body of the orbit, for example the Earth
# * The position and velocity vectors or the orbital elements
# *The time at which the orbit is defined
# 

# First of all, we need to import the relevant modules and classes. Remember, you need to install these modules in your environment before starting running this notebook.

# In[2]:


from astropy import units as u
import numpy as np
from astropy.time import Time 
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit


# # From position and velocity

# There are several methods available to create Orbit objects. For example, if you have the position and velocity vectors you can use from_vectors():

# In[15]:


### Data from Curtis, example 4.3

r = [-6045, -3490, 2500] << u.km
v = [-3.457, 6.618, 2.533] << u.km / u.s

orb = Orbit.from_vectors(Earth, r, v)


# Notice a couple of things:
# 
# * Defining vectorial physical quantities using Astropy units is very easy. The list is automatically converted to a astropy.units.Quantity, which is actually a subclass of NumPy arrays.
# * If you display the orbit you just created, you get a string with the radius of pericenter, radius of apocenter, inclination, reference frame and attractor:

# In[11]:


orb


# * If no time is specified, then a default value is assigned:

# In[12]:


orb.epoch


# In[13]:


orb.epoch.iso


# * The reference frame of the orbit will be one pseudo-intertial frame around the attractor. You can retrieve it using the frame property:

# In[7]:


orb.get_frame()


# # Visualization of the orbit

# In[16]:


orb.plot()


# This plot is made in the so called perifocal frame, which means:
# 
# * you’re visualizing the plane of the orbit itself,
# * the $(x)$ axis points to the pericenter, and
# * the $(y)$ axis is turned $90 \mathrm{^\circ}$ in the direction of the orbit.
# 
# The dotted line represents the osculating orbit: the instantaneous Keplerian orbit at that point. This is relevant in the context of perturbations, when the object shall deviate from its Keplerian orbit.

# # From classical orbital elements

# You can also define an Orbit using a set of six parameters called orbital elements. Although there are several of these element sets, each one with its advantages and drawbacks, right now poliastro supports the classical orbital elements:
# 
# * Semimajor axis $(a)$.
# * Eccentricity $(e)$.
# * Inclination $(i)$.
# * Right ascension of the ascending node $(\Omega)$.
# * Argument of pericenter $(\omega)$.
# * True anomaly $(\nu)$.
# 
# In this case, you’d use the method from_classical():

# In[17]:


### Data for Mars at J2000 from JPL HORIZONS
a = 1.523679 << u.AU
ecc = 0.093315 << u.one
inc = 1.85 << u.deg
raan = 49.562 << u.deg
argp = 286.537 << u.deg
nu = 23.33 << u.deg

orb = Orbit.from_classical(Sun, a, ecc, inc, raan, argp, nu)


# In[18]:


orb.plot()


# Notice that whether you create an Orbit from $(r)$ and $(v)$ or from elements you can access many mathematical properties of the orbit:

# In[19]:


orb.period.to(u.day)


# In[25]:


print(orb.r,"\n",orb.v)


# # Moving forward in time: propagation

# Now that you have defined an orbit, you might be interested in computing how is it going to evolve in the future. In the context of orbital mechanics, this process is known as propagation.
# 
# For example, start by importing an example orbit from the International Space Station:

# In[26]:


from poliastro.examples import iss


# In[27]:


iss


# In[28]:


iss.epoch


# In[29]:


iss.nu.to(u.deg)


# In[30]:


iss.n.to(u.deg / u.min)


# In[31]:


time0=['2013-03-18T12:00:00.000']
times=['2022-10-20T17:02:00']
t=Time(times, format='isot',scale='utc')
t0=Time(time0,format='isot',scale='utc')
DeltaT=t.jd[0]-t0.jd[0]
print(DeltaT)


# Using the propagate() method you can now retrieve the position of the ISS after some time:

# In[32]:


iss_today = iss.propagate(DeltaT << u.day)
iss_today.epoch  # Notice you advanced the epoch!


# In[34]:


from poliastro.twobody.sampling import EpochsArray, TrueAnomalyBounds, EpochBounds
from poliastro.util import time_range

start_date = Time("2022-11-11 05:05", scale="utc")
end_date = Time("2022-11-11 07:05", scale="utc")

# One full revolution
ephem1 = iss.to_ephem()

# Explicit times given
ephem2 = iss.to_ephem(strategy=EpochsArray(epochs=time_range(start_date, end_date)))

# Automatic grid, true anomaly limits
ephem3 = iss.to_ephem(strategy=TrueAnomalyBounds(min_nu=0 << u.deg, max_nu=180 << u.deg))

# Automatic grid, epoch limits
ephem4 = iss.to_ephem(strategy=EpochBounds(min_epoch=start_date, max_epoch=end_date))


# In[ ]:




