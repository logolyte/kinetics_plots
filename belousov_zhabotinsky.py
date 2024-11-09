#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
belousov_zhabotinsky.py

Usage: ./belousov_zhabotinsky.py

Description: Models the oscillations in the Belousov-Zhabotinsky reaction.
Uses the Oregonator model. Rate constants do not necessarily represent the rate
constants of real reactions. The reactions included in the model along with the
rate constants used are:
    1. A + Y -> X + P,      k = 1.28
    2. X + Y -> 2 P,        k = 8*10⁵
    3. A + X -> 2 X + 2 Z,  k = 8
    4. 2X -> A + P,         k = 2*10³
    5. Z + B -> 1/3 Y,      k = 1
The letters represent the following compounds:
    A: BrO₃⁻
    B: CH₂(COOH)₂
    X: HBrO₂
    Y: Br⁻
    Z: [Fe(phen)₃]³⁺ (oxidized ferroin)
    P: HOBr

Dependencies:
    numpy, matplotlib.pyplot, kinetics_plots.py
    
Output: Shows two plots with the concentration of reactants, intermediates and
products over time. One plot shows all compounds included in the model and one
only shows the intermediaries X, Y and Z for which the concentration
oscillates.

Date: 2024-11-08
Author: Love Sundin
"""

from kinetics_plots import Reaction, plot_reactions

# Define reactions
r1 = Reaction({"A": 1, "Y": 1}, {"X": 1, "P": 1}, 1.28)
r2 = Reaction({"X": 1, "Y": 1}, {"P": 2}, 8e5)
r3 = Reaction({"A": 1, "X": 1}, {"X": 2, "Z": 2}, 8)
r4 = Reaction({"X": 2}, {"A": 1, "P": 1}, 2e3)
r5 = Reaction({"Z": 1, "B": 1}, {"Y": 1/3}, 1)

# Plot all compounds
figure = plot_reactions([r1, r2, r3, r4, r5], {"A": 0.06, "Z": 0.00002, "B": 0.06}, time_step=1e-3, time_range=[0,600])
figure.show()
# Plot the intermediates X, Y and Z
figure = plot_reactions([r1, r2, r3, r4, r5], {"A": 0.06, "Z": 0.00002, "B": 0.06}, time_step=1e-3, time_range=[0,600], plot = {"X", "Y", "Z"})
figure.show()