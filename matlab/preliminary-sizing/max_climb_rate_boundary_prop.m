function [W2P] = ...
    max_climb_rate_boundary_prop(W2S, ROC_max, LD_max, C_D_0, K, eta_p,...
    design_density)
%MAX_CLIMB_RATE_BOUNDARY_PROP calculates the maximum power loading for a
% desired maximum climb rate at the given wing loading
% Inputs:
%   W2S:            array containing wing loading values to
%                   calculate admissable power loading for [N/m^2]
%   ROC_max:        desired maximum climb rate [m/s]
%   LD_max:         maximum lift-to-drag ratio [-]
%   C_D_0:          lift-independent drag coefficient [-]
%   K:              induced drag coefficient 1/(pi*AR*e) [-]
%   eta_p:          propeller efficiency [-]
%   design_density: design density [kg/m^3]
% Outputs:
%   W2P:            maximum power loading for desired maximum
%                   climb rate [kg/W]

W2P = 1 /  (ROC_max/eta_p +...
            sqrt(2/design_density/sqrt(3*C_D_0/K).*W2S)*1.155/LD_max/eta_p);
end

