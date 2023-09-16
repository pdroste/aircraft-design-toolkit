function [T2W] = ...
    takeoff_run_boundary_jet(W2S, TOFL, C_D_0_TO, C_L_TO, K,...
    design_density, cf_ground, g)
%TAKEOFF_RUN_BOUNDARY_JET calculates the minimum thrust-to-weight ratio at
% reference conditions for the desired take-off run at design density.
% Inputs:
%   T2W:                array containing wing loading values to
%                       calculate admissable thrust-to-weight ratio [N/m^2]
%   TOFL:               desired take-off field length [m]
%   V_TO:               desired take-off velocity [m/2]
%   C_D_0_TO:           lift-independent drag coefficient
%                       in take-off configuration [-]
%   C_L_TO:             lift coefficient in take-off configuration [-]
%   K:                  induced drag coefficient 1/(pi*AR*e) [-]
%   eta_prop:           propeller efficiency [-]
%   design_density:     air density for design point[kg/m^3]
%   cf_ground:          ground friction coefficient [-]
%   g:                  gravitational acceleration [m/s^2]
% Outputs:
%   T2W:            thrust-to-weight ratio for desired take-off run [kg/W]

% total ground drag coefficient
C_D_g = C_D_0_TO + K * C_L_TO^2 - cf_ground*C_L_TO;

% take-off exponential
exp_TO = exp(0.6 * design_density * g * C_D_g * TOFL ./  W2S);

% admissable thrust to weight
T2W = (cf_ground - (cf_ground + C_D_g/C_L_TO) .* exp_TO) ./ (1 - exp_TO);

end