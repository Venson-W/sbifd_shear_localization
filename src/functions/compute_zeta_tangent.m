function zeta_tangent = compute_zeta_tangent(z, gamma0, tauf, A, sigman, p, eta, dy_trap)
    % Compute tangent of zeta parameter for fault mechanics
    % Inputs:
    %   z           - Depth [m]
    %   gamma0      - Reference strain rate [1/s]
    %   tauf        - Shear stress [Pa]
    %   A           - a-b parameter (dimensionless)
    %   sigman      - Normal stress [Pa]
    %   p           - Pore pressure [Pa]
    %   eta         - Viscosity [Pa⋅s]
    %   yy_gammadot - Strain rate array [1/s]
    % Output:
    %   zeta_tangent - Dimensionless tangent (slope)
    
    % Integrate strain rate derivative over yy_gammadot using trapezoidal rule
    % unit Pa / (m/s) * 1 / (s*Pa) * m
    s = eta*fast_trapz(dstrainratedtau(z, gamma0, tauf, A, sigman, p), dy_trap);
    

    % Compute final result
    zeta_tangent = 1 + s;
end