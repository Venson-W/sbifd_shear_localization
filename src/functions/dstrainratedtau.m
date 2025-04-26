% input format
% strainrate(gamma0(i_dsp,:), tau, sigma_e, a(i_dsp,2:end-1),b(i_dsp,2:end-1), fr0)
% z = fr0 / (a - b)
% A = a - b
function dgammadotdtauf = dstrainratedtau(z, gamma0, tauf, A, sigman, p)
    % Derivative of strain rate with respect to tauf
    % Inputs:
    %   z      - Normalized depth (dimensionless)
    %   gamma0 - Reference strain rate [1/s]
    %   tauf   - Shear stress [Pa]
    %   A      - a-b parameter (dimensionless)
    %   sigman - Normal stress [Pa]
    %   p      - Pore pressure [Pa]
    % Output:
    %   dgammadotdtauf - Derivative of strain rate w.r.t. shear stress [1/(s⋅Pa)]
    
    dgammadotdtauf = 2 .* exp(-z) .* gamma0 .* cosh(tauf ./ (A .* (sigman - p))) ...
        ./ (A .* (sigman - p));
end 

