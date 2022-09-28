#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 12:02:12 2022

@author: konstantinos

Assigmnent 3 and 4 for Tut3
"""
#%% Import AMUSE & declare constants

from amuse.units import units
import amuse.lab as al
import matplotlib.pyplot as plt
import numpy as np # Bound to be useful
from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.community.ph4.interface import ph4

# We're plotting CDFs all the time.
def cdf_plot(data,color=None,label=None,xlabel=None):
    '''
    Parameters
    ----------
    data : Array of floats/ints. Depends on mumpy and matplot
    
    Optionals
    color: str, a matplotlib color
    label: str, a label for the data
    xlabel: str, axis labels
        
    Returns
    -------
    Plots the Empirical CDF of the Data.
    '''
    x  = np.sort(data)
    y  = np.arange(len(data))/float(len(data))
    plt.plot(x,y,
             color=color,
             label=label)
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel('CDF')
    # plt.xscale('log')
    # plt.yscale('log')
    return 0
#%% Make a cluster
amuse_data = 0
def king(n_stars, W0, alpha_IMF):
    # Initial Mass function
    # n_stars = 100
    # alpha_IMF = -2.35 # Exponent of the Salpeter
    smallest = 0.1 | units.MSun # Smallest star in dist
    biggest = 100 | units.MSun # Biggest
    
    m_stars = al.new_powerlaw_mass_distribution(
                        n_stars, # How many stars
                        smallest, # Smallest
                        biggest, # Biggest
                        alpha_IMF # Power law exponent
                        ) 
    # Get this out now, for CDF later. Clunky.
    global amuse_data 
    amuse_data = m_stars.value_in(units.MSun)
    # Tidal radius of the GC
    r_cluster = 1.0 | units.parsec
    
    # Convert to Nbody
    converter = nbody_system.nbody_to_si(m_stars.sum(),r_cluster)
    
    # W0 = Φ(0)/σ^2, thus the bigger it is the stronger the biding.
    # W0 = 3.0
    
    # Make a King model and initialize the converter
    bodies = new_king_model(n_stars, W0, convert_nbody=converter)
    bodies.scale_to_standard(converter)
    return bodies, converter
#%% Evolve it.

# Let's try 2 different w
w_range=[3,6]

# How many times we are going to find the first binary
Repeats = 10

# Lists to put stuff in
bin_list = []
time_list = []

for w in range(len(w_range)):    
    for i in range(Repeats):
        # Let's make this code more robust. While loop
        flag = True
        
        # Make a model
        n_stars = 100
        W0 = w_range[w]
        alpha_IMF = -2.35
        bodies, converter = king(n_stars, W0, alpha_IMF)
        
        # Select gravity solver. ph4 in this case
        gravity = ph4(converter)
        
        # Add our particles to the solver
        gravity.particles.add_particles(bodies)
        
        # Connect our copy with the solver's
        channel = gravity.particles.new_channel_to(bodies)
        
        # Stepping through time
        time = 0 | units.Myr
        step = 0.1 | units.Myr
        
        while flag == True:
            # Do a step
            gravity.evolve_model(time)
            
            # Extract from solver
            channel.copy()
            
            # Time report
            if time.value_in(units.Myr) % 10.0 == 0.0:
                print("Time is: ", time.in_(units.Myr))
                
            # Binary Check
            b = bodies.get_binaries()
            if(len(b)>0):
                flag = False # Stop the loop
                bin_list.append(b)
                time_list.append(time.value_in(units.Myr))
                
            # Time Step
            time = time + step
            
        # Close the channel
        gravity.stop()

            
#%% Make cdf
plt.figure(0)
cdf_plot(time_list[:Repeats],'gold','W:' + str(w_range[0]))
cdf_plot(time_list[Repeats:], 'seagreen','W:' + str(w_range[1]), 'Time of first binary emergence [Myr]')

plt.figure(1)
# The primary and the secondary have the same mass ??
bin_numbers = [ bin_list[i][0][0].mass.value_in(units.MSun) for i in range(len(bin_list))]

# IMF CDFs
cdf_plot(amuse_data, 'red', 'IMF','Binary Mass [MSun]')
# Simulation CDFs
cdf_plot(bin_numbers[:Repeats],'gold','W:' + str(w_range[0]), 'Binary Mass [MSun]')
cdf_plot(bin_numbers[Repeats:], 'seagreen','W:' + str(w_range[1]), 'Binary Mass [MSun]')