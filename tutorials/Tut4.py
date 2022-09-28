
"""
Created on Tue Sep 20 19:33:30 2022

@author: konstantinos

AMUSE Tutorial 4
Nemesis and the Sun
"""

#%% Imports
# Load in the amuse units module, the particle module and 
# some generator for producing some conditions.
from amuse.units import units, constants
from amuse.lab import Particles
from amuse.ext.solarsystem import new_solar_system
from amuse.ext import orbital_elements as oe
import matplotlib.pyplot as plt
import numpy as np

# Make a solar system
system = new_solar_system()

# Give COLORS
colarr = [ 'yellow', 'rosybrown','wheat',
          'royalblue', 'firebrick' , 'orange',
          'khaki', 'cyan' , 'navy' ,'grey' ]
system[0].zorder = 3

for i in range(len(system)):
    system[i].color = colarr[i]
    system.shape = 'o'
    

def plotsystem(solarsystem,title):
    plt.figure()
    for i in range(len(solarsystem)):
        
        if solarsystem[i].zorder == None:
            solarsystem[i].zorder = 1
        
        plt.scatter(solarsystem.x[i].value_in(units.AU),
                    solarsystem.y[i].value_in(units.AU),
                    color = solarsystem.color[i],
                    marker = solarsystem.shape[i],
                    label = solarsystem.name[i],
                    zorder = solarsystem.zorder[i]
                    )
    # plt.legend(loc='best')  
    plt.grid()
    plt.xlabel('x-position [AU]')
    plt.ylabel('y-position [AU]')
    plt.xlim(-1e5, 1e5)
    plt.ylim(-1e5, 1e5)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.title(title)

def plot_ae(primary,secondary):
    orbels = oe.get_orbital_elements_from_arrays(primary.position - secondary.position,
                                                 primary.velocity - secondary.velocity,
                                                 primary.mass+secondary.mass)
    a = orbels[0][0].value_in(units.AU)
    e = orbels[1][0]
    plt.scatter(a,e, 
                color = secondary.color)
    plt.title('Semi Major Axis vs Eccentricity')
    plt.ylabel('Eccentricity')
    plt.xlabel('Semi-Major Axis [AU]')
    plt.ylim(-0.05,1)
    plt.xlim(0,1e5)
    

#%% Add Nemesis

msun = system[0].mass
mnemesis = 0.2 | units.MSun

a = 95000 | units.au
e = 0.7
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
sun_and_nemesis = new_binary_from_orbital_elements(msun, 
                                                    mnemesis, 
                                                    a, e,
                                                    G=constants.G)
sun_and_nemesis[0].name = 'Solar System'
sun_and_nemesis[0].color = 'gold'
sun_and_nemesis[0].shape = 'o'
sun_and_nemesis[1].name = 'Nemesis'
sun_and_nemesis[1].color = 'purple'
sun_and_nemesis[1].shape = 'o'
sun_and_nemesis[1].zorder = 2

# #%% Join the Binary with the Solar System

# # Extract sun from solar system
sun = system[system.name=="SUN"]
# Move the sun to the origin
system.position -= sun.position
system.velocity -= sun.velocity
# Extract solar system and nemesis from binary
ss = sun_and_nemesis[sun_and_nemesis.name=="Solar System"]
nemesis = sun_and_nemesis[sun_and_nemesis.name=="Nemesis"]
# Give to the solar system the attributes of the binary
system.position += ss.position
system.velocity += ss.velocity
# Add Nemesis
system.add_particle(sun_and_nemesis[1])
plotsystem(system, 'Solar system + Nemesis')

# # #%% Assignment 1. Orbit from cartesian


# orb_el = oe.get_orbital_elements_from_arrays(sun.position - nemesis.position,
#                                               sun.velocity - nemesis.velocity,
#                                               sun.mass + nemesis.mass)
# # Question 1. 
# '''
# Used sun's attributes instead of the entire solar
# system for the positions of the primary. Solvable
# problem, just use the proper ones'
# '''
#%% Assignment 2.1 Make the Oort Cloud.
# This is the second implementation
# The first one used generate_binaries
# to create 100 binaries and then system.move_to_center()
# and oort.position += sun.position

# Remove the Terrestial Planets
system.remove_particle( system[1:10] )

# Radii array
oort_radii = np.linspace(1e4, 5e4, 100) | units.AU
oort_mass = 0 | units.kg
oort_cloud = Particles() 

# Make oort cloud
for i in range(len(oort_radii)):
   
    sun_oort_obj_binary = oe.new_binary_from_orbital_elements(
                                        msun, 
                                        oort_mass, 
                                        oort_radii[i], 
                                        eccentricity = 0, # Circular Orbits
                                        G=constants.G)

    # Remove anything the sun might add
    # sun_oort_obj_binary.position -= sun_oort_obj_binary[0].position
    # sun_oort_obj_binary.velocity -= sun_oort_obj_binary[0].velocity

    # Name it oort_obj, for clarity
    oort_obj = sun_oort_obj_binary[1]

    # Add to system
    oort_cloud.add_particle(oort_obj)

# Move to barycenter
oort_cloud.position += system[0].position
oort_cloud.velocity += system[0].velocity

# Add the oort cloud
system.add_particles(oort_cloud)
# system.move_to_center()

for i in range(2,len(system)):
    system[i].color = 'lightsteelblue'
    system[i].shape = '^'
    system[i].name = 'Oort Object: ' + str(i)

# system.move_to_center()
plotsystem(system, 'Check out the Oort cloud') # Plots the correct thing!

# Eccentricity Check
plt.figure()
for i in range(1, len(system)):
    primary = system[0]
    secondary = system[i]
    plot_ae(primary,secondary)
4#%% Assignment 2.2 Set up the solver

# I'll use ph4 because this is the only one i have any experience with
from amuse.community.ph4.interface import ph4

# We need to do it in Nbody units.
from amuse.units import nbody_system
# The biggest mass and the biggest distance should be 1
# so the mass of the sun and the furthest object in the oort cloud
# In python, we are able to index from the end, so [-1] accesses the last object
converter = nbody_system.nbody_to_si(system[0].mass.sum(),system[-1].position.sum())

# Initialize solver
gravity = ph4(converter)
# Add our particles to the solver
# system.scale_to_standard(converter)
gravity.particles.add_particles(system)
# Make a channel
channel = gravity.particles.new_channel_to(system)
nem10 = 160
end = nem10
times = np.arange(0, end, 0.1) | units.Myr
plottimes = np.arange(0, end, 0.2) | units.Myr
#%% Assignment 2.3 Do the thing
j=0
import time # To measure time
import datetime as dt # for now.time
start = time.time()

print('Simulation started')
for i in times:
    # Do a step
    gravity.evolve_model(i)
    # Extract from solver
    channel.copy()
    # Print checks
    print("Sim time: ", i.in_(units.Myr))
    print('Time elapsed: ', time.time()-start, ' seconds')
    now = dt.datetime.now()
    print('Time on print: ', now.time() )

    
    if i in plottimes:
        plotsystem(system, 'At time: ' + str(i.in_(units.Myr)) )
        plt.savefig('graphs/nem2'+str(j)+'.png')
        j=j+1
        # plt.cla()
        
        print('I should print something...')
        
gravity.stop()
#%% Final plot and Ecc check
plotsystem(system, 'Final')

# Semi-major axis - Ecc plot
# Check which has bigger binding energy with sun or Nemesis
# and plot accordingly
plt.figure()
for i in range(len(system)):
    primary = system[0]
    secondary = system[i]
    plot_ae(primary,secondary)
    
#%% Make the .gif
from PIL import Image

i=0
image_paths=[None]*len(plottimes)
im=[]
for i in range(len(image_paths)):
    image_paths[i]='graphs/nem2'+str(i)+'.png' #load up the image names
    new_im= Image.open(image_paths[i]) #open one
    im.append(new_im)

im[0].save('graphs/nem2.gif', format='GIF',
            append_images=im[1:],
            save_all=True,
            duration=500, #in ms
            loop=1) #for zero repeats forever

