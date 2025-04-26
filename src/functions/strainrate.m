% input format
% strainrate(gamma0(i_dsp,:), tau, sigma_e, a(i_dsp,2:end-1),b(i_dsp,2:end-1), fr0)
% z = fr0./(a - b)
% A = a - b
% sigma_e = sigma_zz - p 
function gammadot = strainrate(z, gamma0, tauf, A, sigman, p)
    gammadot = 2 .* exp(-z) .* gamma0 .* sinh(tauf ./ (A .* (sigman - p)));
end 

