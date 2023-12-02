"""
aerodynamics.py

Contains methods related to calculating the aerodynamic properties of a preliminary design.
"""

import numpy as np


def lifting_line(theta, c_l_alpha, chord, total_alpha, halfspan):
    """
    Return the sine Fourier coefficients produced by analysing the lift distribution across a wing using the method
    of lifting line theory.
    :param theta:           polar angles of lift distribution of wing segment (excluding 0 and pi) [rad]
    :param c_l_alpha:       slopes of lift coefficient curve of wing segments [1/rad]
    :param chord:           chord length of wing segments [m]
    :param total_alpha:     total angle of attack of wing segment, including free-stream angle of attack,
                            geometric twist and segment zero-lift angle of attack [rad]
    :param halfspan:        halfspan of the wing [m]
    :return:                uneven Fourier coefficients of the lift distribution [-]
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
    :param fourier_coefficients:    uneven Fourier coefficients of the lift distribution [-]
    :param aspect_ratio:            aspect ratio of the wing [-]
    :return:                        total lift coefficient of the wing [-]
    """
    return np.pi * fourier_coefficients[0] * aspect_ratio


def fourier_induced_drag_coefficient(fourier_coefficients, aspect_ratio):
    """
    Return the total induced drag coefficient of the wing based on the Fourier coefficients of the lift distribution.
    :param fourier_coefficients:    uneven Fourier coefficients of the lift distribution [-]
    :param aspect_ratio:            aspect ratio of the wing [-]
    :return:                        total induced drag coefficient of the wing [-]
    """
    # create uneven coefficient multiplier
    n = np.arange(0, fourier_coefficients.size, 1) * 2 + 1
    n = n.reshape(fourier_coefficients.shape)

    # return induced drag coefficient
    return np.pi * aspect_ratio * np.sum(n * fourier_coefficients**2)


def fourier_span_efficiency_factor(fourier_coefficients):
    """
    Return the span efficiency factor (Oswald factor) of the wing based on the Fourier coefficients of the lift
    distribution.
    :param fourier_coefficients:    uneven Fourier coefficients of the lift distribution [-]
    :return:                        span efficiency factor [-]
    """
    # create uneven coefficient multiplier
    n = np.arange(0, fourier_coefficients.size, 1) * 2 + 1
    n = n.reshape(fourier_coefficients.shape)

    # calculate sum of squared Fourier coefficients
    delta = np.sum(n*fourier_coefficients**2)

    delta /= fourier_coefficients[0]

    # return span efficiency factor
    return 1/(1 + delta)


def fourier_induced_downwash(fourier_coefficients, theta, tas):
    """
    Return the induced downwash distribution of the wing based on the uneven Fourier coefficients of the lift
    distribution.
    :param fourier_coefficients:    uneven Fourier coefficients of the lift distribution [-]
    :param theta:                   polar angles of lift distribution of wing segment [rad]
    :param tas:                     true airspeed [m/s]
    :return:                        downwash distribution of the wing [m/s]
    """
    # create uneven coefficient multiplier
    n = np.arange(0, fourier_coefficients.size, 1) * 2 + 1
    n = n.reshape(fourier_coefficients.shape)

    # calculate induced downwash
    w = n * fourier_coefficients * np.sin(n * theta) / np.sin(theta)

    w *= tas

    # return downwash
    return w


def elliptic_lift_distribution_center(halfspan, lift_desired, rho, tas):
    """
    Returns the product of lift coefficient and chord length at the center of a wing with elliptical lift distribution.
    :param halfspan:        half-span of the wing [m]
    :param lift_desired:    desired total lift of the wing [N]
    :param rho:             air density [kg/m^3]
    :param tas:             true airspeed [m/s]
    :return:                product of lift coefficient and chord length at center of the wing [m]
    """
    return lift_desired/(0.5 * halfspan * rho * tas**2 * np.pi)
