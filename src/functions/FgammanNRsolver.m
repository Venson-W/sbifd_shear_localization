
function gamma = FgammanNRsolver(ny, eta, gamma0, f0, a, b, tau0, sigma_xz, sigma0, sigma_n, pa, p_pre, dy, tol, max_iter)
    % Solve for gamma using Newton-Raphson method
    %
    % Inputs:
    %   ny       - Number of points in y direction (ny = 200)
    %   eta      - Viscosity coefficient
    %   gamma0   - Initial guess for shear rate
    %   f0, a, b - Material constants
    %   tau0     - Stress component
    %   sigma_xz - Shear stress
    %   sigma0, sigma_n, pa, p_pre - Pressure and stress terms
    %   dy       - Grid spacing in y-direction
    %   tol      - Convergence tolerance
    %   max_iter - Maximum iterations
    %
    % Output:
    %   gamma - Solved shear rate values

    % Initialize gamma with an initial guess (small values)
    gamma = gamma0 * ones(ny, 1);
    
    % Constants for computation
    C = 2 * gamma0 / exp(f0 / (a - b));  % Precompute constant
    S_denom = sigma0 + sigma_n - pa - p_pre;  % Denominator for sinh argument

    for iter = 1:max_iter
        % Compute the function F(gamma)
        gamma_sum = sum(gamma(1:end-1) + gamma(2:end));  % Summation term
        S = (tau0 + sigma_xz - eta * (dy / 2) * gamma_sum) / S_denom / (a - b);
        F = gamma - C * sinh(S);
        
        % Check for convergence
        if norm(F, inf) < tol
            disp(['Converged in ', num2str(iter), ' iterations.']);
            return;
        end

        % Construct sparse Jacobian matrix J
        J = speye(ny);  % Identity matrix (diagonal ones)
        coshS = cosh(S);
        coeff = -C * coshS * (-eta * dy / (2 * S_denom));

        % Fill the Jacobian matrix with banded structure
        for i = 1:ny
            J(i, i) = 1 + coeff;  % Main diagonal
            if i < ny
                J(i, i+1) = coeff;  % Upper diagonal
                J(i+1, i) = coeff;  % Lower diagonal
            end
        end

        % Solve for update Δγ using sparse solver
        delta_gamma = -J \ F;

        % Update gamma
        gamma = gamma + delta_gamma;

        % Check for small updates
        if norm(delta_gamma, inf) < tol
            disp(['Newton-Raphson converged in ', num2str(iter), ' iterations.']);
            return;
        end
    end

    disp('Newton-Raphson did not converge within max iterations.');
end
