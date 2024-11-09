#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kinetics_plots.py

Usage: import kinetics_plots

Description: Module for plotting concentrations of compounds over time in
systems of chemical reactions. Contains a Reaction object and a function for
plotting. Uses numpy for calculations and matplotlib.pyplot for plotting.

Dependencies:
    numpy, matplotlib.pyplot
    
Objects:
    Reaction: A class for storing chemical reactions.
    
Functions:
    plot_reactions: Takes a list of Reaction objects and plots the
    concentration of compounds over time.
    a set of 

Date: 2024-11-09
Author: Love Sundin
"""

import numpy
from matplotlib import pyplot

# A class for storing chemical reactions
class Reaction:
    def __init__(self, reactant_dict: dict, product_dict: dict,
                 forward_rate: float, reverse_rate: float = 0):
        """
        Constructor for the Reaction class, which stores a chemical reaction
        for kinetics calculations. Sets the reactants, products, forward rate
        and optionally reverse rate of the reaction.

        Parameters
        ----------
        reactant_dict : dict
            Dictionary with reactants as keys and coefficients as values.
            Reactants are represented by strings.
        product_dict : dict
            Dictionary with products as keys and coefficients as values.
            Products are represented by strings.
        forward_rate : float
            The rate constant for the forward reaction.
        reverse_rate : float, optional
            The rate constant for the reverse reaction. The default is 0.

        Returns
        -------
        None.

        """
        self.reactant_dict = reactant_dict
        self.product_dict = product_dict
        self.forward_rate = forward_rate
        self.reverse_rate = reverse_rate
        
    def __str__(self) -> str:
        """
        Returns a string representation of the reaction.

        Returns
        -------
        return_string : str
            A string showing the reactants and products of the reaction with
            coefficients.

        """
        # The string to return
        return_string = ""
        # Display all reactants
        for reactant_index, (reactant, coefficient) in enumerate(sorted(self.reactant_dict.items())):
            if reactant_index != 0:
                return_string += " + "
            return_string += f"{coefficient} {reactant}"
        return_string += " -> "
        # Display all products
        for product_index, (product, coefficient) in enumerate(sorted(self.product_dict.items())):
            if product_index != 0:
                return_string += " + "
            return_string += f"{coefficient} {product}"
        return return_string
        
    def __repr__(self) -> str:
        """
        Returns a representation of the object that can be used to reconstruct
        it.

        Returns
        -------
        str
            A string showing how the constructor of the object can be called to
            give an identical Reaction object.

        """
        return_string = f"kinetics_plots.Reaction({repr(self.reactant_dict)}, {repr(self.product_dict)}, {self.forward_rate}, {self.reverse_rate})"
        return return_string
        
def plot_reactions(reaction_list: list, initial_concentration_dict: dict,
                   time_step: float = 0.01, time_range: list = [0, 1000],
                   concentration_unit: str = "M", time_unit: str = "s",
                   constant_concentrations: set = set(), plot: set = None,
                   return_matrices: bool = False):
    """
    Plots the concentration of reactants and products over time in a system of
    reactions. Uses Euler's method to integrate concentrations over time.

    Parameters
    ----------
    reaction_list : list
        A list with Reaction objects containing the reactions to plot.
    initial_concentration_dict : dict
        Dictionary where keys are strings corresponding to compounds and values
        are initial concentrations. Compounds not given an initial value have
        their initial value set to 0.
    time_step : float, optional
        The size of each time step. The default is 0.01.
    time_range : list, optional
        The time range to calculate concentrations in. The default is
        [0, 1000].
    concentration_unit : str, optional
        The concentration unit to display. The default is "M".
    time_unit : str, optional
        The time unit to display. The default is "s".
    constant_concentrations : set, optional
        Set of strings corresponding to compounds whose concentration should be
        kept constant. The default is set().
    plot : set, optional
        Set of strings correspoding to compounds to include in the plot. If it
        is None, all reactants and products are plotted. The default is None.
    return_matrices : bool, optional
        Whether or not arrays with concentrations and time at each time step
        should be returned. The default is False.

    Returns
    -------
    pyplot.figure
        A figure where compound concentrations are plotted over time. If
        return_matrices is True, returns a tuple of (figure: pyplot.figure,
        concentration_matrix: numpy.array, time_array: numpy.array) where
        concentration_matrix contains the concentration of each compound after
        each time step and time_array contains the time at each time step. In
        concentration_matrix, rows are time steps and each plotted compound has
        one column. The compounds are given indices in alphabetical order.

    """
    
    # Create a list with all compounds in the reaction
    compound_set = {}
    for reaction in reaction_list:
        compound_set = compound_set | reaction.reactant_dict.keys()
        compound_set = compound_set | reaction.product_dict.keys()
    compound_list = []
    for compound in sorted(compound_set):
        compound_list.append(compound)
    
    compounds = len(compound_list)
    reactions = len(reaction_list)
    
    # List of starting concentrations
    concentration_list = []
    for compound in compound_list:
        if compound in initial_concentration_dict:
            concentration_list.append(initial_concentration_dict[compound])
        else:
            concentration_list.append(0)
    concentration_array = numpy.array(concentration_list, dtype=float)
    
    # Matrix where rows are reactions and columns are the coefficient of
    # compounds in that reaction. Used to calculate the reaction rate at each
    # time point
    reaction_power_matrix = numpy.zeros((2 * reactions, compounds))
    # Matrix where rows are compounds and columns are the reactions,
    # corresponding to how much the reaction affects the concentration of the
    # compound. Used to update the concentration of compounds at each time
    # point
    compound_reaction_matrix = numpy.zeros((compounds, 2 * reactions))
    
    # Set the values of these matrices so they can be used to calculate
    # differentials
    for reaction_index, reaction in enumerate(reaction_list):
        # Consider all forward and reverse reactions
        forward_reaction_index = 2 * reaction_index
        reverse_reaction_index = 2 * reaction_index + 1
        for compound_index, compound in enumerate(compound_list):
            # Iterate over all compounds, and if the compound is a reactant or
            # product update the matrices
            if compound in reaction.reactant_dict:
                # The compound drives the reaction forward
                reaction_power_matrix[forward_reaction_index, compound_index] = reaction.reactant_dict[compound]
                # The compound is consumed in the reaction
                if not compound in constant_concentrations:
                    compound_reaction_matrix[compound_index, forward_reaction_index] += -reaction.forward_rate * reaction.reactant_dict[compound]
                    compound_reaction_matrix[compound_index, reverse_reaction_index] += reaction.reverse_rate * reaction.reactant_dict[compound]
            if compound in reaction.product_dict:
                # The compound drives the reverse reaction
                reaction_power_matrix[reverse_reaction_index, compound_index] = reaction.product_dict[compound]
                # The compound is produced in the reaction
                if not compound in constant_concentrations:
                    compound_reaction_matrix[compound_index, forward_reaction_index] += reaction.forward_rate * reaction.product_dict[compound]
                    compound_reaction_matrix[compound_index, reverse_reaction_index] += -reaction.reverse_rate * reaction.product_dict[compound]
    
    # Make an array of time points
    time_array = numpy.arange(*time_range, time_step)
    
    # Matrix with the concentration of all compounds
    concentration_matrix = []
    
    # Use Euler's method to integrate concentrations over time
    for time in time_array:
        reaction_products = numpy.prod(numpy.power(concentration_array, reaction_power_matrix), axis=1)
        compound_differentials = reaction_products @ compound_reaction_matrix.transpose()
        concentration_array += compound_differentials * time_step
        concentration_matrix.append(concentration_array.copy())
        
    concentration_matrix = numpy.array(concentration_matrix)
    
    # Make a list containing the columns of the compounds to plot
    if not plot is None:
        plot_column_list = []
        for compound in sorted(plot):
            plot_column_list.append(compound_list.index(compound))
        concentration_matrix=concentration_matrix[:, plot_column_list]
    else:
        plot_column_list = list(range(compounds))
    
    # Create a figure and plot concentrations over time
    figure, axis = pyplot.subplots(figsize = (8, 6))
    for column, compound_index in enumerate(plot_column_list):
        compound = compound_list[compound_index]
        axis.plot(time_array, concentration_matrix[:, column], label = compound)
    axis.legend(frameon = False, loc="lower center", bbox_to_anchor=(0.5, 1), ncol=4)
    axis.set_xlabel(f"Time [{time_unit}]")
    axis.set_ylabel(f"Concentration [{concentration_unit}]")
    
    # Return the figure, and if specified by the user also the concentration
    # matrix and time array.
    if return_matrices:
        return figure, concentration_matrix, time_array
    else:
        return figure