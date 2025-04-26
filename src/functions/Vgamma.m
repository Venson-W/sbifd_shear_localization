% The function to evaluate the shear strain rate convergence
function FV = Vgamma(V, z, tauf, sigma_e, a, b, gamma0, dy_trap)

    gammadot = 2 .*gamma0 .* exp(-z) .* sinh(tauf ./ (sigma_e .* (a- b)));
    FV = abs((V - fast_trapz(gammadot, dy_trap) + eps) ./ (V + eps));
    
end