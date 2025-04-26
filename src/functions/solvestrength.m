function [tauf_return, nb_it_return, delta_tauf_return] = solvestrength(tauf, tau_dyn, gamma0, A, sigman, eta, p, z, dy_trap, tol, maxit)
    % Initialize variables
    % Inputs:
    %   tauf        - Initial shear stress [Pa]
    %   tau_dyn     - Dynamic shear stress [Pa]
    %   gamma0      - Reference strain rate [1/s]
    %   A           - a-b parameter (dimensionless)
    %   sigman      - Normal stress [Pa]
    %   eta         - Viscosity [Pa⋅s]
    %   p           - Pore pressure [Pa]
    %   z           - Depth [m]
    %   yy_gammadot - Integration points [1/s]
    %   tol         - Tolerance for convergence
    %   maxit       - Maximum iterations
    
    delta_tauf = 1.1 * tol;
    nb_it = 0;

    % Iterative loop (Newton-Raphson method)
    while (abs(delta_tauf) > tol) && (nb_it < maxit)  
        delta_zeta = -compute_zeta(z, gamma0, tauf, A, sigman, p, eta, dy_trap) + tau_dyn;
        zeta_p = compute_zeta_tangent(z, gamma0, tauf, A, sigman, p, eta, dy_trap);
        delta_tauf = delta_zeta / zeta_p;
        tauf = tauf + delta_tauf;
        nb_it = nb_it + 1;
    end

    % If max iterations reached, display warning
    if nb_it == maxit
        fprintf('Isaac and Joseph failed at finding an equilibrated solution (error %.5f)\n', abs(delta_tauf));
    end
    tauf_return = tauf;
    nb_it_return = nb_it;
    delta_tauf_return = delta_tauf;
end
