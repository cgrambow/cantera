#!/usr/bin/env python
# -*- coding: utf-8 -*-


# Import Mechanism
mechanism = 'chem_annotated.cti'

# make plots?
plot = True
    
# Define initial state of reactor T (K), D (kg/m**3)
T = 438.15 # K
#P = 1.0 # atm
D = 600 #kg/m**3

# names of initial reactants in the cti file
O2_name = 'oxygen(3)'
basestock_name = 'C16H32(4)'

# concentration of initial reactants
O2_c0 = 0.001
basestock_c0 = 1 - O2_c0
p
# Sensitivity analysis
sensitivity_species = ['C16H32(4)']

# Initial time step (seconds)
dt = 1

# write output every d_out seconds
d_out = 100

#write screen output every d_s_out seconds
d_s_out = 10

# number of hours to run
run_time_hours = 24
