"""
aerodynamics.py

Contains methods related to calculating the aerodynamic properties of a preliminary design.
"""

import numpy as np


def lifting_line(theta, c_l_alpha, chord, total_alpha, halfspan):
    """
    Returns the sine Fourier coefficients produced by analysing the lift distribution across a wing using the method
    of lifting line theory.
    :param theta:           polar angles of lift distribution of wing segment
    :param c_l_alpha:       slopes of lift coefficient curve of wing segments
    :param chord:           chord length of wing segments
    :param total_alpha:     total angle of attack of wing segment, including free-stream angle of attack,
                            geometric twist and segment zero-lift angle of attack
    :param halfspan:        halfspan of the wing
    :return:                uneven Fourier coefficients of the lift distribution
    """

    # uneven multiplicator
    n = 2 * np.arange(0, theta.size, 1).reshape((1, theta.size)) + 1

    # construct Fourier matrix
    fourier = np.sin(theta) * np.ones((1, theta.size))

    fourier += c_l_alpha * chord/8/halfspan * n

    fourier *= np.sin(theta * n)

    # construct circulation matrix
    circulation = c_l_alpha * chord/8/halfspan * np.sin(theta) * total_alpha

    return np.linalg.solve(fourier, circulation)


def fourier_total_lift_coefficient(fourier_coefficients, aspect_ratio):
    """
    Return the total lift coefficient of the wing based on the first uneven Fourier coefficient.
    :param fourier_coefficients:    uneven Fourier coefficients of the lift distribution
    :param aspect_ratio:            aspect ratio of the wing
    :return:                        total lift coefficient of the wing
    """
    return np.pi * fourier_coefficients[0] * aspect_ratio

def fourier_induced_drag_coefficient(fourier_coefficients, aspect_ratio):
    """
    Return the total induced drag coefficient of the wing based on the Fourier coefficients of the lift distribution.
    :param fourier_coefficients:    uneven Fourier coefficients of the lift distribution
    :param aspect_ratio:            aspect ratio of the wing
    :return:                        total induced drag coefficient of the wing
    """
    #create uneven coefficient multiplier
    n = np.arange(0, fourier_coefficients.size, 1) * 2 + 1

    # return induced drag coefficient
    return np.pi * aspect_ratio * np.sum(n * fourier_coefficients**2)

