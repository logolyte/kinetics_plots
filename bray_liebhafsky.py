#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bray_liebhafsky.py

Usage: ./bray_liebhafsky.py

Description: Models the oscillations in the Bray-Liebhafsky reaction. The model
used is based on the model described in "A mathematical model of the
Bray–Liebhafsky reaction" by Dimsey et al. from 2024
(https://doi.org/10.1098/rspa.2023.0964). Rate constants used do not
necessarily represent the rate constants of real reactions. The reactions
included in the model along with the rate constants used are:
    1. A + Y -> X + P,      k = 3.5*10⁻³
    2. X + Y -> 2 P,        k = 10
    3. B + X -> 2 X + Z,    k = 4*10⁻³
    4. 2 X -> P + A,        k = 0.1
    5. Z -> Y,              k = 1*10⁻²
    6. Z -> Q,              k = 1*10⁻⁴
The letters represent the following compounds:
    A: IO₃⁻
    B: H₂O₂
    X: HIO₂
    Y: I⁻
    Z: I₂ (aq)
    P: HOI
    Q: I₂ (g)

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
r1 = Reaction({"A": 1, "Y": 1}, {"X": 1, "P": 1}, 3.5e-3)
r2 = Reaction({"X": 1, "Y": 1}, {"P": 2}, 10)
r3 = Reaction({"B": 1, "X": 1}, {"X": 2, "Z": 1}, 4e-3)
r4 = Reaction({"X": 2}, {"P": 1, "A": 1}, 1e-1)
r5 = Reaction({"Z": 1}, {"Y": 1}, 1e-2)
r6 = Reaction({"Z": 1}, {"Q": 1}, 1e-4)

# Plot all compounds
figure = plot_reactions([r1, r2, r3, r4, r5, r6], {"A": 10, "B": 10, "X": 0.001, "Y": 0.001, "Z": 0.001}, time_range = [0, 1500])
figure.show()
# Plot only X, Y and Z
figure = plot_reactions([r1, r2, r3, r4, r5, r6], {"A": 10, "B": 10, "X": 0.001, "Y": 0.001, "Z": 0.001}, plot={"X", "Y", "Z"}, time_range = [0, 1500])
figure.show()