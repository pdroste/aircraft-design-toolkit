function [W2S] = ...
    stall_speed_boundary(V_stall, C_L_max, design_density)
%STALL_SPEED_BOUNDARY calculates the maximum wing loading for a desired
% stall speed.
% Inputs:
%   V_stall:            desired stall speed [m/s]
%   C_L_max:            maximum lift coefficient [-]
%   design_density:     air density at design point[kg/m^3]
% Outputs:
%   W2S:                wing loading [N/m^2]

W2S = design_density/2 * V_stall.^2 * C_L_max;

end

