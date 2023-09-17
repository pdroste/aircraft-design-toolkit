function [T2W] = ...
    max_speed_boundary_jet(W2S, V_max, ref_density, design_density,...
                        C_D0, K)
%MAX_SPEED_BOUNDARY_JET calculates the minimum thrust-to-weight ratio at
% reference conditions for the desired maximum speed at design density.
% Inputs:
%   T2W:            array containing wing loading values to
%                   calculate admissable thrust-to-weight ratio [N/m^2] for
%   V_max:          desired maximum speed [m/s]
%   ref_density:    reference density [kg/m^3]
%   design_density: design density [kg/m^3]
%   C_D0:           lift-independent drag coefficient [-]
%   K:              induced drag coefficient 1/(pi*AR*e) [-]
% Outputs:
%   T2W:            thrust-to-weight ratio for desired maximum speed [kg/W]

density_ratio = design_density/ref_density;

a = 0.5 * ref_density * C_D0;
b = 2 * K / design_density / density_ratio;

T2W = a * V_max.^2 ./ W2S + b ./ V_max.^2 .* W2S;

end

