"""
aerodynamics.py

Contains methods related to calculating the aerodynamic properties of a preliminary design.
"""

import numpy as np


def lifting_line(theta, c_l_alpha, chord, total_alpha, halfspan):

    # uneven multiplicator
    n = 2 * np.arange(0, theta.size, 1).reshape((1, theta.size)) + 1

    # construct Fourier matrix
    fourier = np.sin(theta) * np.ones((1, theta.size))

    fourier += c_l_alpha * chord/8/halfspan * n

    fourier *= np.sin(theta * n)

    # construct circulation matrix
    circulation = c_l_alpha * chord/8/halfspan * np.sin(theta) * total_alpha

    return np.linalg.solve(fourier, circulation)
