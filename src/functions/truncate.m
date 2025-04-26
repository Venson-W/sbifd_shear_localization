function [M_kernel, ts_index] = truncate(M_kernelTV, nx, nt_max, c, k, dt, t_tol, v, vu, kernel_name, ODE_fitter_switch)
    M_kernel = cell(1, nx);
    ts_index = zeros(1,length(k));
    if  ODE_fitter_switch == 0
        if kernel_name == "diffusion" % K2
            for i_xx = 1:nx
               for i_tt = 1:nt_max
                    % M_kernelT = -c * k(i_xx).^2 * exp(-c * k(i_xx) ^ 2 * (i_tt*dt - dt/2)) / (sqrt(pi * (i_tt*dt - dt/2) * c * k(i_xx).^2));
                    M_kernelT = -exp(-c .* k(i_xx) ^ 2 * (i_tt*dt - dt/2)) / sqrt((i_tt*dt - dt/2));
                    % longer allocation time
                    if i_tt == 1       
                        M_kernelTV(i_tt,1) = M_kernelT;
                        ts_index(1,i_xx) = i_tt;
                    % if the tolerance reached, update the truncating point and store the updated kernel 
                    elseif M_kernelT / M_kernelTV(1) > t_tol
                        ts_index(1,i_xx) = i_tt;
                        M_kernelTV(i_tt,1) = M_kernelT;    
                    else
                        ts_index(1,i_xx) = i_tt;
                        M_kernel{i_xx}(i_tt,1) = M_kernelT; 
                       continue 
                    end
               end
               M_kernel{i_xx}(:,1) = M_kernelTV(1:i_tt,1);
            end
        else if kernel_name == "delta_x" %K1
           for i_xx = 1:nx
               for i_tt = 1:nt_max
     
                    M_kernelT = -erfc(sqrt(c * k(i_xx).^2 * (i_tt*dt - dt/2)));
                    % longer allocation time
                    if i_tt == 1       
                        M_kernelTV(i_tt,1) = M_kernelT;
                        ts_index(1,i_xx) = i_tt;
                    % if the tolerance reached, update the truncating point and store the updated kernel 
                    elseif M_kernelT / M_kernelTV(1) > t_tol
                        ts_index(1,i_xx) = i_tt;
                        M_kernelTV(i_tt,1) = M_kernelT;    
                    else
                        ts_index(1,i_xx) = i_tt;
                        M_kernel{i_xx}(i_tt,1) = M_kernelT; 
                       continue 
                    end
               end
               M_kernel{i_xx}(:,1) = M_kernelTV(1:i_tt,1);
            end
        end
        end
    end
   
end