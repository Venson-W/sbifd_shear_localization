function [F_gamma_result, p]  = Fgammapsolver(sigma_xz, tau0, G, vs, yy_gamma, gamma, Tn, Te1d, dt, pn, pe1d, gamman, ...
    sigma_zz0, sigma_zz, pa, gamma0, fr0, a, b, i_dsp, j, i_t)
    
    % Recalculate tau, T, and p based on the current gamma
    V_dsp = trapz(yy_gamma, gamma, 2);
    tau = sigma_xz(i_dsp,j) + tau0(i_dsp,j) - G/(2*vs)*repmat(V_dsp, 1, length(yy_gamma));

    % impose tau(tau<0) = 0;
    %if any(tau(:)<0)
    %    % Compute coefficient where tau < 0
    %    negative_tau_mask = tau < 0;  % Identify where tau is negative
   % 
   %     denominator = (G/(2*vs)) * repmat(V_dsp, 1, length(yy_gamma));
   %     penalty_coef = negative_tau_mask.*((sigma_xz(i_dsp,j) + tau0(i_dsp,j)) ./ denominator);
   %     penalty_coef(penalty_coef == 0) = 1;  % Ensure coef does not exceed 1
        
        % Adjust gamma
    %    gamma = gamma .* penalty_coef;
        
        % correct to make min(tau,[],'all') = 0
    %    tau = sigma_xz(i_dsp,j) + tau0(i_dsp,j) - G/(2*vs)*repmat(V_dsp, 1, length(yy_gamma));
    %end 

    T(i_dsp,j) = Tn(i_dsp,j) + Te1d(i_t*dt, Tn, gamma, tau, i_dsp) * dt;
    delta_T = (T(i_dsp,j) - Tn(i_dsp,j)) / dt;
    p(i_dsp,j) = pn(i_dsp,j) + pe1d(i_t*dt, pn, gamman, gamma, delta_T, dt, i_dsp) * dt;
    
    
    
    sigma_e = sigma_zz0 + sigma_zz(i_dsp,j) - p(i_dsp,j) - pa;
    
    %sigma_e(sigma_e < 0) = 1e6;
    
    % Calculate F(gamma)
    F_gamma_result = gamma - 2 .*gamma0(i_dsp,:) .* sinh(tau ./ (sigma_e .* (a(i_dsp,j) - b(i_dsp,j)))) ...
                     ./ (exp(fr0./(a(i_dsp,j) - b(i_dsp,j))));
    
end