#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 13:18:25 2022

@author: konstantinos
AMUSE Tutorial #3: Simple N-body

"""
#%% Import AMUSE & declare constants

from amuse.units import units
import amuse.lab as al
import matplotlib.pyplot as plt
import numpy as np # Bound to be useful

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
    plt.xscale('log')
    # plt.yscale('log')
    return 0
#%% Testing the Salpeter IMF
n_stars = 100
alpha_IMF = -2.35 # Exponent of the Salpeter
smallest = 0.1 | units.MSun # Smallest star in dist
biggest = 100 | units.MSun # Biggest

m_stars = al.new_powerlaw_mass_distribution(
                    n_stars, # How many stars
                    smallest, # Smallest
                    biggest, # Biggest
                    alpha_IMF # Power law exponent
                    ) 

# CDF for the generated stars by AMUSE
amuse_data = m_stars.value_in(units.MSun)
cdf_plot(amuse_data, 'red', 'amuse')

# Define a theoretical Salpeter dist.
def salpeter(x):
    sal = np.zeros(len(x)-1)
    for i in range(len(x)-1):
        dm = x[i+1] - x[i]
        sal[i] = x[i]**(-2.35) * dm
    return sal

mass_range = np.linspace(0.1,100,1000)
theoretical_data = salpeter(mass_range)
cdf_plot(theoretical_data,'blue','Theory','Stellar Masses [M_sol]')

#%% Converter

r_cluster = 1.0 | units.parsec
from amuse.units import nbody_system
converter = nbody_system.nbody_to_si(m_stars.sum(),r_cluster)
from amuse.ic.kingmodel import new_king_model
# W0 = Φ(0)/σ^2, thus the bigger it is the stronger the biding.
W0 = 3.0
bodies = new_king_model(n_stars, W0, convert_nbody=converter)
bodies.scale_to_standard(converter)

#%% Check what the converter did
from amuse.plot import scatter # This apparently has built in cmap support

def plot_snapshot(bodies):
    v = (bodies.vx**2 + bodies.vy**2 + bodies.vz**2).sqrt()
    scatter(bodies.x, bodies.y, c=v.value_in(units.kms), alpha=0.5)
    plt.colorbar()
    
plot_snapshot(bodies)

#%% Actually doing it

from amuse.community.ph4.interface import ph4
from amuse.ext.LagrangianRadii import LagrangianRadii

# Select gravity solver. ph4 in this case
gravity = ph4(converter)
# Add our particles to the solver
gravity.particles.add_particles(bodies)

# A channel is a 'permanent' connection to a code's particle
# set. Multiple calls to a code's particle set need to set up
# a new connection every time; with a channel, we can copy
# information back without opening a new connection.
# This does not automatically update bodies! See below
channel = gravity.particles.new_channel_to(bodies)

times = np.arange(0, 100, 0.1) | units.Myr
RL25 = [] | units.parsec
Rvir = [] | units.parsec

for time in times:
    # Do a step
    gravity.evolve_model(time)
    # Extract from solver
    channel.copy()
    # Lists for Virial, Lagrange and RL25(?)
    Rvir.append(bodies.virial_radius())
    L = LagrangianRadii(bodies)
    RL25.append(LagrangianRadii(bodies)[5])

    if not time.value_in(units.Myr)%10.0:
        print("cluster at Time=", time.in_(units.Myr), 
              "Mass=", bodies.mass.sum().in_(units.MSun),
              "Rvir=", Rvir[-1].in_(units.parsec))
    b = bodies.get_binaries()
    if(len(b)>0):
        print("Number of binaries found:", len(b))
        print("first binary:", b[0])


plt.plot(times.value_in(units.Myr), RL25.value_in(units.parsec))
plt.plot(times.value_in(units.Myr), Rvir.value_in(units.parsec))
plt.show()

#%% Assignment 2. w8 for binary

channel = gravity.particles.new_channel_to(bodies)

# Just make time huge and add a break
times = np.arange(0, 1e5, 0.1) | units.Myr
RL25 = [] | units.parsec
Rvir = [] | units.parsec

# For the cumulative read a bit on king models, change it
# add another for loop to store the data and run a bunch of times
# say 30. No it is not statisticallty significant.

# Also store the masses of the binaries

for time in times:
    # Do a step
    gravity.evolve_model(time)
    # Extract from solver
    channel.copy()
    # Lists for Virial, Lagrange and RL25(?)
    Rvir.append(bodies.virial_radius())
    L = LagrangianRadii(bodies)
    RL25.append(LagrangianRadii(bodies)[5])

    if not time.value_in(units.Myr)%10.0:
        print("cluster at Time=", time.in_(units.Myr), 
              "Mass=", bodies.mass.sum().in_(units.MSun),
              "Rvir=", Rvir[-1].in_(units.parsec))
    b = bodies.get_binaries()
    if(len(b)>0):
        print("Number of binaries found:", len(b))
        print("first binary:", b[0])
        


plt.plot(times.value_in(units.Myr), RL25.value_in(units.parsec))
plt.plot(times.value_in(units.Myr), Rvir.value_in(units.parsec))
plt.show()


