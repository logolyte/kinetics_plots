#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lotka_voltera.py

Usage: ./lotka_voltera.py

Description: Models the oscillations in a chemical reaction based on
Lotka-Voltera models of predator-prey interactions. The reaction mechanism is
completely theoretical and does not represent the mechanism of any known
reaction. The reactions included in the model along with the rate constants
used are:
    1. A + X -> 2 X, k = 0.06
    2. X + Y -> 2 Y, k = 0.6
    3. Y -> B,       k = 0.06

Dependencies:
    numpy, matplotlib.pyplot, kinetics_plots.py
    
Output: Shows a plot with the concentration of all compounds over time.

Date: 2024-11-08
Author: Love Sundin
"""

from kinetics_plots import Reaction, plot_reactions

# Define reactions
r1 = Reaction({"A": 1, "X": 1}, {"X": 2}, 0.06)
r2 = Reaction({"X": 1, "Y": 1}, {"Y": 2}, 0.6)
r3 = Reaction({"Y": 1}, {"B": 1}, 0.06)
# Plot the change in concentrations over time
figure = plot_reactions([r1, r2, r3], {"A": 8, "X": 0.1, "Y": 0.05})
figure.show()