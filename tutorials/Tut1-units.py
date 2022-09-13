# -*- coding: utf-8 -*-
"""
AMUSE Tutorial #1: Scalar Units
"""

#%% Import AMUSE & declare constants

from amuse.units import units
from amuse.units.constants import G
#%% Assignment 1
# Earth Orbital Velocity

# The Earth orbits the Sun, usually.
mstar = 1 | units.MSun # Mass is 1 Solar Mass

# At a distance of 1 AU, approximately.
distance = 1 | units.AU

# The Equation is v=sqrt( G M / d )
v_orb = (G*mstar/distance).sqrt()

# Print, please
print( 'Orbital Velocity of the Earth:' , v_orb.in_(units.kms))

#%% Assignment 2
# S2 Escapes from Sgr A*

# Sgr A* is a massive ball
m_BH = 4.154e6 | units.MSun 

# S2 Pericenter Distance
pericenter = 120 | units.AU

# Escape Velocity at the Pericenter
v_esc = (2*G*m_BH/pericenter).sqrt()

# Print, please
print( 'Escape Velocity of Sgr A* at S2\'s pericenter: \n' , v_esc.in_(units.kms))

#%% Question 1
# Terror! Asteroids impact the Earth.

import numpy as np # Need it for ranges

# It does orbit the sun

mstar = 1 | units.MSun

# If the Semi-Major is smaller than 0.5 AU
# then the Asteroid will not intersect with
# Earths orbit. Thus, we are safe.


# Ellipticity range

min_el = 0
max_el = 0.4
el_range = np.linspace(min_el, max_el, 50)


r_earth = 1 | units.AU

"""
The minumum orbital velocity will be when the collison occurs at the aphelion
where r=1 AU = a(1-e)
So the aphelion velocity reduces to
"""
v_min = np.sqrt(G*mstar* (1-max_el)/r_earth )

"""
The minumum orbital velocity will be when the collison occurs at the perihelion
where r=1 AU = a(1+e)
So the perihelion velocity reduces to:
"""
v_max = np.sqrt(G*mstar* (1+max_el)/r_earth )
        
print('The range of orbital velocities of meteors that can impact the earth is: \n', 
      v_min.value_in(units.kms), '[Km/s] -', v_max.value_in(units.kms), '[Km/s]')

'''
 This happens after L1. So from 0.01 AU the asteroid falls toward the earth'
'''

mplanet = 1 | units.MEarth 
rplanet = 1 | units.REarth 
fall_distance = 0.01 | units.AU 
v_ff = np.sqrt(2*G*mplanet * (1/rplanet - 1/(fall_distance+rplanet) ))

print('The gravitational field of the Earth provides an additional: \n', 
       v_ff.value_in(units.kms), '[Km/s]')

'''
Asumming the orbital velocity does not change during the free fall
from verctror addition we get the following
'''

v_max = np.sqrt(v_max**2 + v_ff**2)
v_min = np.sqrt(v_min**2 + v_ff**2)

print('Thus the range of velocities is: \n',
      (v_min).value_in(units.kms), '[Km/s] -', (v_max).value_in(units.kms), '[Km/s]')

#%% Question 2
# Lux Confusion

from amuse.units.constants import Stefan_hyphen_Boltzmann_constant

sigma = Stefan_hyphen_Boltzmann_constant
T_wiki = 5772 | units.K
F_wiki = sigma * T_wiki**4

Sun_radius = 1 | units.RSun
L_wiki = F_wiki * 4 * np.pi * Sun_radius**2
L_amuse = 1 | units.LSun

print('Wikipedia Luminosity',(L_wiki).value_in(units.LSun), '[LSun] \n'
       'Amuse Luminosity', (L_amuse).value_in(units.LSun), '[LSun]')

F_amuse = L_amuse / (4*np.pi*Sun_radius**2)
T_amuse = (F_amuse / sigma)**(1/4)

print('Wikipedia Temp.',(T_wiki).value_in(units.K), '[LSun] \n'
       'Amuse Temp', (T_amuse).value_in(units.K), '[LSun]')


