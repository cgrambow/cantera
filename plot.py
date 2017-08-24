# -*- coding: utf-8 -*-

from rmgpy.tools.plot import GenericPlot, SimulationPlot

def plot_profiles(data, topSpecies=10):
    """
    data is tuple:
    (time data, [temperature data, density data, species1 data, species2 data, ...])
    """

    time, data_list = data
    temperature = data_list[0]
    density = data_list[1]
    species_data = [d for d in data_list[2:] if d.index != -1]  # remove unwanted species

    # plot
    GenericPlot(xVar=time, yVar=temperature).plot('temperature.png')
    GenericPlot(xVar=time, yVar=density).plot('density.png')
    SimulationPlot(xVar=time, yVar=species_data, numSpecies=topSpecies, ylabel='Mole Fraction').plot('mole_fractions.png')
