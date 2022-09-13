#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 12:31:08 2022
@author: konstantinos

AMUSE Tutorial #2: Particles

"""
#%% Import AMUSE & declare constants

from amuse.units import units
from amuse.lab import Particles
import numpy as np # Bound to be useful
#%% Making the Sun and the Earth

# Make two particles
sun_and_earth = Particles(2)
# First one is the sun
sun = sun_and_earth[0]
# It has m,x,v
sun.mass = 1 | units.MSun 
sun.position = (0,0,0) | units.au 
sun.velocity = (0,0,0) | units.kms 
# print("Sun=", sun)

# Second one is the earth
earth = sun_and_earth[1]
# It has m,x,
earth.mass = 1 | units.MEarth
earth.position = (1,0,0) | units.au

# We need to calculate velocity
from amuse.units.constants import G

def rel_orb_vel(m,d):
    v = np.sqrt(G*m/d)
    return v

v_orb_earth = rel_orb_vel( sun_and_earth.mass.sum(), earth.position.sum())
earth.velocity = (0,1,0) * v_orb_earth
# print("Earth=", earth)

# Move them to the center of Mass
sun_and_earth.move_to_center()

#%% Attributes

# Let's name them, inside amuse
setattr(sun_and_earth, "name", "")
sun_and_earth.name = ['Sun', 'Earth']

#%% The moooon

moon = Particles(1)
moon.name = "moon"
moon.mass = 7.34767309e+22 | units.kg
moon.position = (384400, 0, 0) | units.km
vorb = rel_orb_vel(earth.mass + moon.mass, 
                                 moon.position.sum())
moon.velocity = (0, 1, 0) * vorb
# print('Moon= ',moon)

# So that's the moon. We need it to be around the earth though, not around
# the sun which we have placed at the origin

moon.position = moon.position + earth.position
moon.velocity = moon.velocity + earth.velocity

# Now we can place it in our system
sun_and_earth.add_particle(moon)

# Perhaps also rename it
solarsystem = sun_and_earth
solarsystem.move_to_center()

#%%  Assignment 1: Zeus

jupiter = Particles(1)
jupiter.name = 'Jupiter'
jupiter.mass = 1 | units.MJupiter 
jupiter.position = (5.2038,0,0) | units.AU # Assuming a circular orbit
vorb = rel_orb_vel(jupiter.mass + sun.mass, jupiter.position.sum())
jupiter.velocity = (0,1,0) * vorb

solarsystem.add_particle(jupiter)

#%% Assignment 2: Inclinations & Anomalies

from amuse.ext import orbital_elements

# Let us generate random numbers for the inclination.
# a = 0 to b = Pi/6. We multiply by (b - a) * random() + a
np.random.seed(8)
rng = np.random.default_rng()
inclinations = np.pi/6 * rng.random(2) | units.rad # Realistic inclinations

mean_anomalies = 2*np.pi * rng.random(2) | units.rad # Just anywhere along the orbit

# Mean -> Eccentric
# Assuming circular orbits, the Mean Anomaly IS the True Anomaly
# Use M = E - e * sinE, if needed

# Eccentric -> True. Assume Circular Orbis.
true_anomalies = orbital_elements.true_anomaly_from_eccentric_anomaly(mean_anomalies, 0)

# This gets us the positions and velocities from masses-postions-inclinations
# true anomalies.

posvel = orbital_elements.rel_posvel_arrays_from_orbital_elements(
        primary_mass=sun.mass,
        secondary_mass=earth.mass,
        semi_major_axis=earth.position.sum(),
        eccentricity=0,
        true_anomaly=true_anomalies[0],
        inclination=inclinations[0],
        )

earth.position = posvel[0][0]
earth.velocity = posvel[1][0]

# Do agian for moon
# Remember we are only looking at the Earth-Moon system so reinitialize
# position and velocity

moon.position = (384400, 0, 0) | units.km # Redefine
vorb = rel_orb_vel(earth.mass + moon.mass, 
                                 moon.position.sum())
moon.velocity = (0, 1, 0) * vorb

posvel = orbital_elements.rel_posvel_arrays_from_orbital_elements(
        earth.mass,
        moon.mass,
        moon.position.sum(),
        0,
        true_anomalies[1],
        inclinations[1],
        G=G
        )

moon.position = posvel[0][0] + earth.position # Moon orbits earth
moon.velocity = posvel[1][0] + earth.velocity 

#%% Assignment 3. Binding energy

def binding_energy(m1,m2,a):
    E = G*m1*m2/(2*a) 
    return E


def be_system(solarsystem):
    """
    Parameters
    ----------
    solarsystem : Particle Set. A system of objects whose total
    binding energy we are calculating.

    Returns
    -------
    be, float the binding energy

    """
    be_sum = 0 | units.J
    for i in range(len(solarsystem)):
        for j in range(len(solarsystem)):
            if j<=i:
                continue
            
            be_sum +=  binding_energy(solarsystem.mass[i],solarsystem.mass[j], 
                            solarsystem.distances_squared(solarsystem[i]).sqrt()[j])

    return be_sum

print('Binding Energy:', be_system(solarsystem) ,'[J]')    
# No matter where we move it, it will not change    
solarsystem.position += 100 | units.parsec
solarsystem.vz += 100 | units.kms  
print('Binding Energy:', be_system(solarsystem) ,'[J]') 
solarsystem.move_to_center()
#%% Question 1. Get Binaries

binaries = solarsystem.get_binaries()
print(binaries)
# Binaries goes weird if I move the solar system. Maybe it has to do
# with the hardness???

#%% Assignment 4
system2 = Particles(3)

# New Star
star = system2[0]
star.mass = 2 | units.MSun 
star.position = (0,0,0) | units.AU
star.velocity = (0,0,0) | units.kms

# Planet 1
planet1 = system2[1]
planet1.mass = 10 | units.MEarth
planet1.position = (0.1,0,0) | units.AU
vorb = rel_orb_vel( star.mass + planet1.mass, planet1.position.sum())
planet1.velocity =  (0,1,0) * vorb 

# Planet 2
planet2 = system2[1]
planet2.mass = 100 | units.MEarth
planet1.position = (0.6,0,0) | units.AU
vorb = rel_orb_vel( star.mass + planet2.mass, planet2.position.sum())
planet2.velocity = vorb * (0,1,0) 

# Add it
newsystem = orbital_elements.new_binary_from_orbital_elements(solarsystem.mass.sum(),
                                             system2.mass.sum(),
                                             60 |units.AU,
                                             eccentricity=0.6|units.rad,
                                             true_anomaly=np.pi |units.rad)

newsystem.move_to_center()
#%% Question 2

# I do not understand the question.
# Are there many orbits? I can calculate the planet with the highest binding energy.



















