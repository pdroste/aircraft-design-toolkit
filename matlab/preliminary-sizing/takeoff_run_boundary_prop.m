function [W2P] = ...
    takeoff_run_boundary_prop(W2S, TOFL, V_TO, C_D_0_TO, C_L_TO, K,...
    eta_prop, design_density, cf_ground, g)
%TAKEOFF_RUN_BOUNDARY_PROP calculates the maximum admissable power loading
%for the given wing loading
% Inputs:
%   W2S:                array containing wing loading values to
%                       calculate admissable power loading for [N/m^2]
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
%   W2P:                maximum power loading for desired take-off
%                       characteristics

% total ground drag coefficient
C_D_g = C_D_0_TO + K * C_L_TO^2 - cf_ground*C_L_TO;

% take-off exponential
exp_TO = exp(eta_prop * design_density * g * C_D_g * TOFL ./  W2S);

% admissable power loading
W2P = 1 - exp_TO ./ (cf_ground - (cf_ground + C_D_g/C_L_TO) .* exp_TO) .*...
    eta_prop ./ V_TO;

end

