#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description of what this does.
"""


# Standard library
import cPickle as pickle
import numpy as np
import re

# Dependencies
import cantera as ct
from rmgpy.species import Species
from rmgpy.tools.data import GenericData

# Modules
import batchReactor_functions as cthf
from plot import plot_profiles

# User input
from user_input import *


def initialize_reactor(gas, sensitivity_species):
    num_reactions = gas.n_reactions
    r = ct.IdealGasReactor(gas, energy='off')
    sim = ct.ReactorNet([r])

    if sensitivity_species:
        for i in xrange(num_reactions):
            r.add_sensitivity_reaction(i)
        sim.rtol_sensitivity = 1e-4
        sim.atol_sensitivity = 1e-6

    return r, sim

##############################################################################
##############################################################################
##############################################################################

# import the gas model and set the initial conditions
gas = ct.Solution(mechanism)

# set Temp, Density, Composition of initial mixture
gas.TDX = T, D, 'C16H32(4):{0},oxygen(3):{1}'.format(basestock_c0,O2_c0)

# get the names of all the species
species_names = gas.species_names
# create a look to find the index of oxygen so we can reset it,
# thus implimenting constant concentration of oxygen
n_mech = len(species_names)
index_O2 = [i for i,s in enumerate(species_names) if ('oxygen(3)') == s]

# Intial mass fraction for O2 and N2
O2_0 = gas.X[index_O2]

################################################################################
print('Setting state of reactor')

# Ideal Gas Reactor
# r = ct.IdealGasReactor(gas,energy='off')
# r_vol = r.volume

# sim = ct.ReactorNet([r])
r, sim = initialize_reactor(gas, sensitivity_species)

# Initial time
time = 0.0


# number of seconds to run
run_time = run_time_hours * 3600

# No of steps
#N = run_time / dt
################################################################################

# Array for storing data
#times = np.zeros(N)
#data = np.zeros((N,7))

################################################################################

TDX = gas.TDX

###############   OUTPUT INITIALIZATION #####################################

###################### X vs time output and descriptors
# write the column names as smiles
out_lbl = 'SMILES'
out_dat_SMILES = cthf.label_to_smile(gas.species_names)
out_write = np.append(out_lbl,out_dat_SMILES)
cthf.output_writer(out_write,'./out_data_species_X.csv')


# now add a row of exact MW
out_lbl = 'exact_mass'
out_dat = cthf.rd_kit_exact_mass(out_dat_SMILES)
out_write = np.append(out_lbl,out_dat)
cthf.output_writer(out_write,'./out_data_species_X.csv')


# now add a row of exact radical electrons
out_lbl = 'num_rad_electrons'
out_dat = cthf.rd_kit_num_rad_electrons(out_dat_SMILES)
out_write = np.append(out_lbl,out_dat)
cthf.output_writer(out_write,'./out_data_species_X.csv')

# Initialize lists for storing condition and species data
out_times = []
out_temperature = []
out_density = []
out_species_data = []

###################### thermo vs time output #######
# write the column names
out_lbl = ['T', 'P', 'rho', 'Xbasestock', 'Xoxygen', 'dt']
out_lbl = np.insert(out_lbl, 0, 'time', 0)
cthf.output_writer(out_lbl,'./out_data_thermo.csv')


###################### elements vs time #########
# write the column names
out_lbl = ['C', 'H', 'O']
out_lbl = np.insert(out_lbl, 0, 'time', 0)
cthf.output_writer(out_lbl,'./out_data_elemental_mol.csv')

###################### elements vs time #########
# write the column names
out_lbl = ['C', 'H', 'O']
out_lbl = np.insert(out_lbl, 0, 'time', 0)
cthf.output_writer(out_lbl,'./out_data_elemental_mass.csv')

while time <= run_time:
    
    # Const Concent: reset the mole fraction of O2 and initial mole fraction
    TDX[2][index_O2] = 0.0
    s = sum(TDX[2]) / (1.0 - O2_0)
    TDX[2][index_O2] = O2_0 * s
    gas.TDX = TDX[0], TDX[1], TDX[2]
    
    # Initialize the reactor    
    # r = ct.IdealGasReactor(gas, energy='off')
    # sim = ct.ReactorNet([r])
    r, sim = initialize_reactor(gas, sensitivity_species)
    
    






    ############## output things ##############################################
    
    ## bad hack since modulo is oversensitive to, e.g., 0.00000000001
    ## and the cumulative effect of using dt leads to drift in time.
    if (np.mod(np.around(time,decimals=5),np.around(d_out,decimals=5)) == 0
        or time==0):
        
        ####### species concentration output
        out_dat = r.thermo.X
        # insert a column for time at column index 0
        out_dat = np.insert(out_dat, 0, time, 0)
        cthf.output_writer(out_dat,'./out_data_species_X.csv')

        ####### store data in memory
        out_times.append(time)
        out_temperature.append(TDX[0])
        out_density.append(TDX[1])
        out_species_data.append(r.thermo.X)

        ######## elemental composition output
        # first copy the gas and remove the O2 species so we only look at O
        # incorporated into molecules
        gas_no_O2 = gas
        TDX_no_O2 = r.thermo.TDX
        TDX_no_O2[2][index_O2] = 0
        gas_no_O2.TDX = TDX_no_O2[0], TDX_no_O2[1], TDX_no_O2[2]
        # get the data on elments to output
        out_dat = [gas_no_O2.elemental_mole_fraction('C'), gas_no_O2.elemental_mole_fraction('H'), gas_no_O2.elemental_mole_fraction('O')]
        # insert a column for time at colum index 0
        out_dat = np.insert(out_dat, 0, time, 0)
        cthf.output_writer(out_dat,'./out_data_elemental_mol.csv')

        ######## elemental composition output
        # first copy the gas and remove the O2 species so we only look at O
        # incorporated into molecules
        ### this was done above

        # get the data on elments to output
        out_dat = [gas_no_O2.elemental_mass_fraction('C'), gas_no_O2.elemental_mass_fraction('H'), gas_no_O2.elemental_mass_fraction('O')]
        # insert a column for time at colum index 0
        out_dat = np.insert(out_dat, 0, time, 0)
        cthf.output_writer(out_dat,'./out_data_elemental_mass.csv')

        ######## thermo output
        out_dat = [time, r.T, r.thermo.P, 1/r.thermo.volume_mass, r.thermo[basestock_name].X[0], r.thermo[O2_name].X[0], dt]
        cthf.output_writer(out_dat,'./out_data_thermo.csv')


    # Advance simulation
    sim.advance(dt)

    
    
    ## bad hack since modulo is oversensitive to, e.g., 0.00000000001
    ## and the cumulative effect of using dt leads to drift in time.
    if (np.mod(np.around(time,decimals=5),np.around(d_s_out,decimals=5)) == 0
        or time==0):    
        # grab some data and output to stdout
        print('%10.5f %10.2f %10.3f %10.5f %10.5f %10.5f' % (time, r.T, 1/r.thermo.volume_mass, r.thermo[basestock_name].X, ((r.thermo[O2_name].X - O2_0)*100 / O2_0), dt))
    
    
    # advance time and store TDX for reinit of gas
    time +=  dt
    TDX = r.thermo.TDX    
    
    
    # this is a corrector that reduced dt if oxygen concentration varried
    # too much between resets
    if (abs((r.thermo[O2_name].X - O2_0)/ O2_0) > 0.1):
        dt = 0.1 * dt
    
    # here we can increase dt of O2 content is relatively constant
    # but only do it at same modulo of output so we are sure that it
    # occurs on an mod=0 timestep so we don't loose output
    if ((abs((r.thermo[O2_name].X - O2_0)/ O2_0) < 0.01)
      and dt < d_out
      and np.mod(np.around(time,decimals=5),np.around(d_out,decimals=5)) == 0):
        dt = 10 * dt

# Convert data
out_species_data = np.array(out_species_data)
time = GenericData(label='Time', data=out_times, units='s')
out_temperature = GenericData(label='Temperature', data=out_temperature, units='K')
out_density = GenericData(label='Density', data=out_density, units='kg/m3')
TDX_data = [out_temperature, out_density]

for i, species_name in enumerate(species_names):
    match = re.search('\([0-9]*\)', species_name)
    if match:
        index = int(match.group()[1:-1])
    else:
        index = -1  # don't plot these species (inerts and non-reactive solvent)

    species_generic_data = GenericData(label=out_dat_SMILES[i],
                                       species=Species().fromSMILES(out_dat_SMILES[i]),
                                       data=out_species_data[:,i],
                                       index=index,
                                    )
    TDX_data.append(species_generic_data)

# Write data to disk
all_data = (time, TDX_data)
with open('all_data.dat', 'w') as f:
    pickle.dump(all_data, f)

# Plot data
if plot:
    plot_profiles(all_data)
