#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Input and settings.
"""


## Input

# Chemkin mechanism
mechanism = 'chem_annotated.inp'

# Define reactor state T (K), D (kg/m3)
T = 438.15
D = 600

# Set run time (hr)
run_time_hours = 24

# Set initial species from mechanism file
O2_name = 'oxygen(3)'
basestock_name = 'C16H32(4)'

# Set initial mole fractions
O2_c0 = 0.001
basestock_c0 = 1 - O2_c0

# Set sensitivity analysis species (or empty list)
sensitivity_species = ['C16H32(4)']

## Settings

# Set plot flag
plot = True

# Set initial time step (s)
dt = 1

# Write output every d_out seconds
d_out = 100

# Write screen output every d_s_out seconds
d_s_out = 10
