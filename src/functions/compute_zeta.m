function zeta = compute_zeta(z, gamma0, tauf, A, sigman, p, eta, dy_trap)
    % Compute zeta parameter for fault mechanics
    % Inputs:
    %   z           - Depth [m]
    %   gamma0      - Reference strain rate [1/s]
    %   tauf        - Shear stress [Pa]
    %   A           - a-b parameter (dimensionless)
    %   sigman      - Normal stress [Pa]
    %   p           - Pore pressure [Pa]
    %   eta         - Viscosity [Pa⋅s]
    %   yy_gammadot - Integration points for numerical quadrature [1/s]
    % Output:
    %   zeta        - Stress [Pa]
    
    % Iterate over w using a loop
    % unit Pa / (m/s) * m/s
    s = eta*fast_trapz(strainrate(z, gamma0, tauf, A, sigman, p), dy_trap);
    
    % Compute final result
    zeta = tauf + s;
end