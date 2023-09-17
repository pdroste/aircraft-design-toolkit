function [T2W] = service_ceiling_boundary_jet(W2S, ROC_ceil, LD_max, ...
    C_D_0, K, ref_density, design_density)
%SERVICE_CEILING_BOUNDARY_PROP calculates the minimum thrust-to-weight for a
% desired ceiling and residual ceiling climb rate at the given wing loading
% Inputs:
%   W2S:            array containing wing loading values to
%                   calculate admissable power loading for [N/m^2]
%   ROC_ceil:       desired residual climb rate at ceiling [m/s]
%   LD_max:         maximum lift-to-drag ratio [-]
%   C_D_0:          lift-independent drag coefficient [-]
%   K:              induced drag coefficient 1/(pi*AR*e) [-]
%   eta_p:          propeller efficiency [-]
%   ref_density:    reference density [kg/m^3]
%   design_density: design density at ceiling [kg/m^3]
% Outputs:
%   T2W:            minimum thrust-to-weight ratio for desired maximum
%                   climb rate [kg/W]


density_ratio = design_density/ref_density;

T2W = ROC_ceil ./ ...
    (density_ratio .* sqrt(2./design_density./sqrt(C_D_0/K).*W2S)) + ...
    1 / density_ratio ./ LD_max;

end
