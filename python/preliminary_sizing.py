"""
preliminary_sizing.py

Contains methods for determining the boundaries of allowable thrust-to-weight ratios (for jet aircraft) or power
loading values (for prop-driven aircraft).
"""

import numpy as np


def stall_speed_boundary(stall_speed, c_l_max, design_density):
    """
    Retuns the maximum wing loading for a desired stall speed at given maximum lift coefficient and design density.
    :param stall_speed:         desired stall speed [m/s]
    :type stall_speed:          float
    :param c_l_max:             maximum lift coeffiicent [-]
    :type c_l_max:              float
    :param design_density:      density at which to achieve stall speed [kg/m^3]
    :type design_density:       float
    :return:                    maximum wing loading [N/m^2]
    :rtype:                     float
    """

    return design_density/2 * stall_speed**2 * c_l_max


def service_ceiling_boundary_jet(wing_loading, roc_ceiling, ld_max, c_d_0, K, reference_density, design_density):
    """
    Returns the minimum thrust-to-weight ratio for a given wing loading array for a desired ceiling and residual
    ceiling climb rate
    :param wing_loading:        wing loading [N/m^2]
    :type wing_loading:         np.ndarray
    :param roc_ceiling:         residual rate of climb at ceiling [m/s]
    :type roc_ceiling:          float
    :param ld_max:              maximum lift-to-drag ratio [-]
    :type ld_max:               float
    :param c_d_0:               lift-independent drag coefficient
    :type c_d_0:                float
    :param K:                   induced drag coefficient 1/pi*AR*e [-]
    :type K:                    float
    :param reference_density:   thrust reference density [kg/m^3]
    :type reference_density:    float
    :param design_density:      design density [kg/m^3]
    :type design_density:       float
    :return:                    minimum thrust-to-weight ratio for given wing loading [-]
    :rtype:                     np.ndarray
    """

    density_ratio = design_density/reference_density

    return roc_ceiling / (density_ratio * np.sqrt(2/design_density/np.sqrt(c_d_0/K)*wing_loading)) + \
        1/density_ratio/ld_max


def service_ceiling_boundary_prop(W2S, roc_ceiling, ld_max,
    c_d_0, K, prop_efficiency, reference_density, design_density):
    """

    :param W2S:                 wing loading [N/m^2]
    :type W2S:                  np.ndarray
    :param roc_ceiling:         desired residual climb rate at ceiling [m/s]
    :type roc_ceiling:          float
    :param ld_max:              maximum lift-to-drag ratio [-]
    :type ld_max:               float
    :param K:                   induced drag coefficient 1/(pi*AR*e) [-]
    :type K:                    float
    :param prop_efficiency:     propeller efficiency [-]
    :type prop_efficiency:      float
    :param reference_density:   reference density [kg/m^3]
    :type reference_density:    float
    :param design_density:      design density [kg/m^3]
    :type design_density:       float
    :return:                    maximum power loading for desired wing loading values [N/W]
    :rtype:                     np.ndarray
    """

    density_ratio = design_density / reference_density

    return density_ratio / \
        (roc_ceiling / prop_efficiency +
            np.sqrt(2 / (design_density * np.sqrt(3 * c_d_0 / K)) * W2S) *
            1.155 / ld_max / prop_efficiency)
