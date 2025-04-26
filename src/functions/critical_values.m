% The function to calculate the characteristic values
% delta_xc: characteristic displacement
% gamma_c: characteristic stress
% V_c: characteristic velocity
% tau_c: characteristic time
% L_c: characteristic length
function [delta_xc, gamma_c, V_c, tau_c, L_c] = critical_values(rho_c, Lambda, h, fr0, G, vs, sigma_zz0, pa, tau_b, fr, alpha_hy, alpha_th)
    
    delta_xc = rho_c*h/(fr0*Lambda);

    gamma_c = 2*vs*(fr0.*(sigma_zz0 - pa))/(h*G);
    
    V_c = 2 * vs * mean(tau_b(:)) / G;

    tau_c = gamma_c*h*G/(2*vs)/1e6;

    L_c = 4 ./ fr.^2 .* (rho_c/Lambda).^2 .* (sqrt(alpha_hy) + sqrt(alpha_th)).^2 ./ V_c;
end
