% The function to evaluate the shear strain rate convergence
function Fgamma = Fgamma(gamma, gamma0, tau, sigma_e, fr0, a, b)
    
    Fgamma = gamma - 2 .*gamma0 .* sinh(tau ./ (sigma_e .* (a(:,2:end-1) - b(:,2:end-1)))) ./ (exp(fr0./(a(:,2:end-1) - b(:,2:end-1))));
    
end