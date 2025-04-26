function velocity = integratevelocity(z, gamma0, tauf, A, sigman, p, dy_trap)
    % Compute slip velocity by integrating strain rate
    % Inputs:
    %   z           - Depth [m]
    %   gamma0      - Reference strain rate [1/s]
    %   tauf        - Shear stress [Pa]
    %   A           - a-b parameter (dimensionless)
    %   sigman      - Normal stress [Pa]
    %   p           - Pore pressure [Pa]
    %   yy_gammadot - Integration points for numerical quadrature [1/s]
    % Output:
    %   velocity    - Slip velocity [m/s]
    
    %velocity = trapz(yy_gammadot, strainrate(z, gamma0, tauf, A, sigman, p), 2);
    velocity = fast_trapz(strainrate(z, gamma0, tauf, A, sigman, p), dy_trap);
    %strainrate(z, gamma0, tauf, A, sigman, p)
end
