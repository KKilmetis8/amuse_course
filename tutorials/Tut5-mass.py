#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 18:41:26 2022

@author: konstantinos
AMUSE Tutorial 5
Stellar Evo, the Notebook
different masses
"""

#%% Imports

from amuse.units import units, constants
from amuse.datamodel import Particles
import matplotlib.pyplot as plt
import numpy as np
from amuse.lab import new_kroupa_mass_distribution
from amuse.lab import new_salpeter_mass_distribution
from amuse.community.seba.interface import SeBa

#%% IMF

n_stars = 1024
mmin = 0.1 | units.MSun
mmax = 100 | units.MSun

# Kroupa Dist
mkroupa = new_kroupa_mass_distribution(n_stars,
                                           mass_min=mmin, 
                                           mass_max=mmax)
k_stars = Particles(mass=mkroupa)
# Salpeter Dist
msalpeter = new_salpeter_mass_distribution(n_stars, 
                                           mass_min=mmin, 
                                           mass_max=mmax)
s_stars = Particles(mass=msalpeter)

#%% Set up the code

def start_stellar_code(stars):
    stellar = SeBa()
    stellar.particles.add_particles(stars)
    channels = {"to_stars": stellar.particles.new_channel_to(stars), 
                "to_stellar": stars.new_channel_to(stellar.particles)}
    return stellar, channels
kstellar, kchannels = start_stellar_code(k_stars)
sstellar, schannels = start_stellar_code(s_stars)

#%% Run it

times = 10**np.arange(0.0, 5.0, 0.1) | units.Myr
mmean = []
for time in times:
    kstellar.evolve_model(time)
    kchannels["to_stars"].copy()
    sstellar.evolve_model(time)
    schannels["to_stars"].copy()
    mmean.append(np.mean(k_stars.mass)/np.mean(s_stars.mass))
kstellar.stop()
sstellar.stop()

#%% Count types

def CO_counter(stars, title):
    ''' Counts compact objects in a particle set, plos a bar plot and returns
    an the result.
    ~
    Input: Particle, set containing stars.
           title, str. The title of the plot
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
                                        
counter1 = CO_counter(k_stars, 'Kroupa')
plt.figure()
counter2 = CO_counter(s_stars, 'Salpeter')    

#%% Plot
   
from amuse.plot import plot
plot(times, mmean)
plt.ylabel("relative mean mass")
plt.semilogx()
plt.show()

from amuse.plot import scatter
scatter(s_stars.temperature, s_stars.luminosity, c='r')
scatter(k_stars.temperature, k_stars.luminosity, c='b', s=3)
plt.xlim(2.e+4, 2000)
plt.ylim(1.e-5, 1000)
plt.loglog()
plt.show()




