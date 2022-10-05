#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 09:51:40 2022

@author: konstantinos
AMUSE Tutorial 5
Stellar Evo
different Z
---
How to add metallicity ?
As a method in the particle | No.
As SeBa parameter after adding particles | No.
As SeBa parameter BEFORE adding particles | Yes!        
"""

#%% Imports

from amuse.units import units
from amuse.datamodel import Particles
import matplotlib.pyplot as plt
import numpy as np
from amuse.lab import new_kroupa_mass_distribution
from amuse.community.seba.interface import SeBa

#%% IMF

n_stars = 1024
mmin = 0.1 | units.MSun
mmax = 100 | units.MSun

# Kroupa Dist
mkroupa = new_kroupa_mass_distribution(n_stars,
                                           mass_min=mmin, 
                                           mass_max=mmax)
z_sol_stars = Particles(mass=mkroupa)
z_sol_tenth_stars = Particles(mass=mkroupa)
#%% Initialize SeBa

def start_stellar_code(stars, metallicity):
    stellar = SeBa()
    stellar.parameters.metallicity = metallicity
    stellar.particles.add_particles(stars)
    channels = {"to_stars": stellar.particles.new_channel_to(stars), 
                "to_stellar": stars.new_channel_to(stellar.particles)}
    return stellar, channels

solar_z = 0.02 # Solar Metalicity

# Call the function and then add metallicity
Z_sol_stellar, Z_sol_channels = start_stellar_code(z_sol_stars, solar_z)

# Call the function and then add metallicity, again.
Z_sol_tenth_stellar, Z_sol_tenth_channels = start_stellar_code(z_sol_tenth_stars, solar_z / 10)

#%% Evolve!

times = 10**np.arange(0.0, 5.5, 0.1) | units.Myr
for time in times:
    # Evolve and extract
    Z_sol_stellar.evolve_model(time)
    Z_sol_channels["to_stars"].copy()
    print('High Z', time)
    # For the other pop as well
    Z_sol_tenth_stellar.evolve_model(time)
    Z_sol_tenth_channels["to_stars"].copy()
    print('Low Z', time)

# Stop the code
Z_sol_stellar.stop()
Z_sol_tenth_stellar.stop()

#%% Count types

def CO_counter(stars, title):
    ''' Counts compact objects in a particle set, plos a bar plot and returns
    an the result.
    ~
    Input: Particle set containing stars.
    ---
    Output: List of ints. Contains the counts of each compact object
    the order is, He WD, C-O WD, O-Ne WD, NS, BH
    '''
    # Compact objects are stored in AMUSE as
    he_wd = 10
    co_wd = 11
    on_wd = 12
    ns = 13
    bh = 14

    counter = np.zeros(5)
    for i in range(len(stars)):
        if stars.stellar_type[i].value == he_wd:
            counter[0] += 1
        if stars.stellar_type[i].value == co_wd:
            counter[1] += 1
        if stars.stellar_type[i].value == on_wd:
            counter[2] += 1
        if stars.stellar_type[i].value == ns:
            counter[3] += 1
        if stars.stellar_type[i].value == bh:
            counter[4] += 1
            
    col_arr = ['navajowhite','palegreen','pink', 'deepskyblue' , 'k']
    labels = ['Helium WDs', 'Carbon-Oxygen WDs', 'Oxygen-Neon WDs', 'Neutron Stars <3',
              'Black Hole']
    f = plt.figure()
    f.set_figwidth(10)
    f.set_figheight(10)
    plt.bar(x = labels, height = counter, color = col_arr)
    plt.title(title)
    return counter
                                        
counter1 = CO_counter(z_sol_stars, 'Solar Z')
plt.figure()
counter2 = CO_counter(z_sol_tenth_stars, '1/10 of Solar Z')            
            
        
#%% Plot

from amuse.plot import scatter
scatter(z_sol_stars.temperature, z_sol_stars.luminosity,
        c='gold', label='Z = 0.2')
scatter(z_sol_tenth_stars.temperature, z_sol_tenth_stars.luminosity,
        c='brown', s=3, label='Z = 0.02')
plt.legend()
plt.xlim(2.e+4, 2000)
plt.ylim(1.e-5, 1000)
plt.loglog()
plt.show()





















