% This code was created 
% Author: Yu-Han W.
% Date: 04 02 2025
% Modified:
% Use relative error instead of Fgamma
% Modified:
% use gamman for Tn and tau in gamma updating
% Modified:
% Changed the Fgammasolver and set gamma(isnan(gamma)) == a limiter value
% Modified:
% Include inf and NaN into abnormal conditions
clear;  
clc; 
close all; % clean history
currentFolder = pwd;
[parentFolder, ~, ~] = fileparts(currentFolder);
outputFolder = mfilename;
if ~exist(fullfile(parentFolder, 'post', outputFolder), 'dir')
    % If the folder does not exist, create it.
    mkdir(fullfile(parentFolder, 'post', outputFolder));
    fprintf('Created folder: %s\n', outputFolder);
else
    fprintf('Folder already exists: %s\n', outputFolder);
    delete(fullfile(parentFolder, 'post', outputFolder, '*'));
    fprintf('Folder has been clear: %s\n', outputFolder);
end
addpath([pwd '/functions']);
addpath([pwd '/postprocessing']);
addpath([pwd '/params']);
freq_post = 4000;
rng(41231232); 
%% Model configuration
inputparser_case1;

% Geometry and grid
L  = parser.L; % Fault length and x-axis grid counts
nx = parser.nx; % Counts of the x-axis node
xx = linspace(-L/2 , L/2, nx); % X-axis node coordination
dx = L/(nx - 1); % X-axis grid size
ny = parser.ny; % Counts of the y-axis node
h  = parser.h; % Model thickness
W_fz = h; % Width of gouge
W_fc = 0.000; % Width of fault core
Boundary_S = L/4 - W_fz/2; % The upper boundary 
Boundary_N = L/4 + W_fz/2; % The lower boundary 
dy = W_fz/(ny-1); % Y-axis grid size
yy = linspace(-h/2, h/2, ny);

% Poroelastic boundary 
% Spectral boundary integral
vu = parser.vu; % Undrained Poisson's ratio
v = parser.v; % Drained Poisson's ratio
G = parser.G; % Shear modulus
alphaBiot = parser.alphaBiot;
B = parser.B;
cs = parser.cs; % Shear wave velocity
vs = cs;
% Gouge parameters from a depth-averaged damage material, Table 1, Platt et al., 2014
% Continuum domain, finite difference method here
Lambda         = parser.Lambda;%parser.Lambda;
rho_c          = parser.rho_c;
alpha_th       = parser.alpha_th;
alpha_hx       = parser.alpha_hx;
alpha_hy       = parser.alpha_hy;
eps_pl         = 0;
fr0            = parser.fr0;
fr             = parser.fr;
beta           = parser.beta; % Compressibility (1/Pa) 
beta_T         = parser.beta_T; % 3e-6
a = parser.a*ones(nx,ny); % Rate-dependent
b = parser.b*ones(nx,ny); % Rate-dependent
% Define permeability field with mobility \kappa
% The value of M_kappa_x greatly affects the rupture behavior.
% We defibe a scaling coef (dx/dy)^2 that scales the charac. diff. time to be
% the same in X and Y directions.
% 1) M_kappa_x = kappa*(dx/dy)^2/1e1, the rupture tip appears to be diffusive
% 1) M_kappa_x = kappa*(dx/dy)^2/1e5, the rupture tip appears to be sharp

c              = alpha_hy; % Diffusivity (m2/s)
kappa          = c * beta; % Mobility (m2/(Pa*s))
M_kappa_x      = kappa * ones(nx+1, ny+1)*(dx/dy)^2/10; % X-direction mobility matrix //1e5    
M_kappa_y      = kappa * ones(nx+1, ny+1); % Y-direction mobility matrix
alpha_hostrock = parser.alpha_hostrock;
kappa_hostrock = alpha_hostrock * beta; % Mobility for the host rock

% Define two layers near the boundary to be the host rock
M_kappa_y(:, 1:2) = kappa_hostrock; % (m2/(Pa*s))
M_kappa_y(:, end-1:end) = kappa_hostrock; % (m2/(Pa*s))
M_kappa_x(:, 1:2) = kappa_hostrock; % (m2/(Pa*s))*(dx/dy)^2/1e2
M_kappa_x(:, end-1:end) = kappa_hostrock; % (m2/(Pa*s))*(dx/dy)^2/1e2
% Define diffusivity matrices 
M_alpha_hx = M_kappa_x/beta; % (m2/s)
M_alpha_hy = M_kappa_y/beta; % (m2/s)

% Define heat conductivity field 
kappa_T        = alpha_th * beta_T;
M_kappa_T_x    = kappa_T * ones(nx+1, ny+1); 
M_kappa_T_y    = kappa_T * ones(nx+1, ny+1); %*(dx/dy)^2/1e2
% Define diffusivity matrices 
M_alpha_T_x    = M_kappa_T_x/beta_T; 
M_alpha_T_y    = M_kappa_T_y/beta_T;

% Initial condition
pa             = 50e6; % Background pore pressure
sigma_zz0      = 126e6; % Background normal stress
sigma_xz       = zeros(nx, ny); % Normal stress due to poroelasticity 
sigma_zz       = zeros(nx, ny); % Normal stress due to poroelasticity
T0             = 0 + 483.15; % Background temperature (K)
normalMatrix   = randn(nx, ny)*1;
% Scale the values to range [0, 1]
minVal         = min(normalMatrix(:)); % Minimum value in the matrix
maxVal         = max(normalMatrix(:)); % Maximum value in the matrix
scaledMatrix   = (normalMatrix - minVal) / (maxVal - minVal);
T              = zeros(nx,ny); % Initial temperature (K)
V0             = 1e-9; % Initial shearing
gamma0         = repmat(ones(1,nx)*V0/h,ny,1)'; % Initial strain rate [s-1]
const_gamma0   = V0/h;
delta_x        = zeros(1,nx); % Initial slip (m)
p              = (sigma_zz0 - pa)*scaledMatrix*0.01*0; %initial p matrix *pa*0.001

p(:,1:2)       = 0;
p(:,end-1:end) = 0;
eta = G/(2*cs);
taub           =  fr0.*(sigma_zz0 - pa);

% Characteristic values
[delta_xc, gamma_c, V_c, tau_c, L_c] = critical_values(rho_c, Lambda, h, fr0, G, cs, sigma_zz0, pa, taub, fr, alpha_hy, alpha_th);
delta_xc
gamma_c
V_c
tau_c
L_c
fr_c = tau_c./taub;
normal_dist = @(xx, xi, pii) pii * exp(-xx.^2/(2*xi).^2);
%p_trigger = (sigma_zz0 - pa)*normal_dist(xx,0.05*L,0.57)';
p = p + (sigma_zz0 - pa)*repmat(normal_dist(xx,0.005*L,0.82)', 1, length(yy)); 


taub0 = ones(nx,1)*fr0.*(sigma_zz0 - pa);
taub = fr0.*(sigma_zz0 - pa); %+ ones(nx,1)* V0 * eta + tau_c*1e6*normal_dist(xx,0.005*L,0.6)'
% Time step determination
i_t            = 0; % Initial time step count
sigma          = parser.sigma; % Relaxation coefficient

t_end_showup   = parser.endingtime * vs*fr0*sigma_zz0/h/G;
t_end_showup
t_charc_diffusion_x = dx^2 / max(M_alpha_hx(:));
t_charc_diffusion_x
t_charc_diffusion_y = dy^2 / max(M_alpha_hy(:));
t_charc_diffusion_y
t_charc_heat_x = dx^2 / max(M_alpha_T_x(:));
t_charc_heat_x
t_charc_heat_y = dy^2 / max(M_alpha_T_y(:));
t_charc_heat_y
dt             = sigma * min([dx^2 / max(M_alpha_hx(:)) ...
                     dy^2 / max(M_alpha_hy(:)) ...
                     min(dx, dy)^2 / alpha_th]); % Time step size (s)
t_end          = parser.endingtime; % Simulation time (s) parser.endingtime parser.endingtime
dt0            = dt; % Time step for SBI domain
nt             = round(t_end / dt);  
cumutime       = 0;
disp(['The total time step count is ', num2str(nt)]); % Unit: 1

% Guarantee homogeneous strain rate across the fault
lambda_shr = 2*pi*sqrt((a - b)*rho_c*(alpha_hy+alpha_th) ...
             ./fr0./Lambda./(fr0 + 2*(a-b))./gamma0);
h_0        = max(lambda_shr,[],'all')/ 2;
if h > h_0
    disp('h is too large for creating uniform strain rate');
    disp(['h_0 = ', num2str(h_0)]);
else
    disp('h value is fine');
end
kernelfac = 1;

% Coupling interface distance to the boundary
% compute k
P = xx(end) - xx(1);
k_x = fftshift(2 * pi * (-nx/2:1:(nx/2 - 1)) ./ P); % the length for k should be 2^n
k_x(1) = k_x(2)/10;
%de_gibbs = sinc((-nx/2:1:nx/2-1)/(nx/2));
filter = (sin(pi*(abs([-nx/2:1:nx/2-1]))/(nx*0.5))./(pi*(abs([-nx/2:1:nx/2-1]))/(nx*0.5))).^(1);
filter(isnan(filter)) = 1;
filter_twist = filter(1:nx/2);
filter(1:nx/2)=filter(nx/2+1:end);
filter(nx/2+1:end)=filter_twist;
%filter = ones(nx,1);
% Pre-compute constant outside the loop to save time 
% p - delta x
akdhat1 = (1 + vu)./(1 - vu) .* (1i .* k_x .* G .* B)./3.*filter;%
% p - delta x convolution
akdhat2 = (1 + vu)/(1 - vu) .* (1i .* abs(k_x) .* k_x .* G .* B)./3 .* filter .* sqrt(kappa_hostrock/beta)./ sqrt(pi);%
% p - fluid flux  convolution
akdhat3 = sqrt(alpha_hostrock) ./ sqrt(pi) ./ kappa_hostrock;
% stress xz - delta x 
akdhat4 = -G .* abs(k_x) / (2*(1 - vu));
% stress xz - delta x convolution
akdhat5 = -G * abs(k_x).^3 / (2*(1 - vu)) * 2*(vu - v)/(1-v) * kappa_hostrock/beta;
% stress xz - fluid flux convolution
akdhat6 = 1i .* 3 .* k_x / (4*B*kappa_hostrock*(1+vu)).*filter * 2*(vu - v)/(1-v) * kappa_hostrock/beta; 
% stress zz - fluid flux  convolution
akdhat7 = 3 .* abs(k_x) ./ (4*B*kappa_hostrock*(1+vu)) * 2*(vu - v)/(1-v) * kappa_hostrock/beta;
% T - heat flux  convolution
akdhatT = sqrt(alpha_th)./ sqrt(pi) ./ kappa_T;   

%% Preallocation
% If set maxn to 1 there shouldn't be any downsampling
% generate the kernel as that of the original code
nt_max = round(t_end/ dt);
t_tol = 1e-3;
M_kernelTV = zeros(nt_max,1);
tol_up = 1.0e-4;
int_tol = 1.0e-4; 
maxn = 15; 

% 4978, t_tol = 1e-3
% Define the folder to save pre-allocated kernels
kernel_folder = strcat('kernel_conv_rupture_', outputFolder);
kernel_name = 'kernel_list.mat';

% Check if the folder exists
if (~isfolder(kernel_folder)) || nt_max > 10e4
    
    % If the folder does not exist, create it
    mkdir(kernel_folder);
    fprintf('Folder "%s" created.\n', kernel_folder);

    % Create some sample data to save in the .mat file
    [M_kernel, ts_index] = truncate(M_kernelTV, nx, nt_max, kappa_hostrock/beta, k_x, dt, t_tol, v, vu, "diffusion", 0);
    [M_kernel_T, ts_index_T] = truncate(M_kernelTV, nx, nt_max, alpha_th, k_x, dt, t_tol, v, vu, "diffusion", 0);
    [M_kernel_deltax, ts_index_deltax] = truncate(M_kernelTV, nx, nt_max, kappa_hostrock/beta, k_x, dt, t_tol, v, vu, "delta_x", 0);
    
    % redue the size of kernel matrix
    % Diffusion kernel - pore pressure, p
    [M_kernelR, indr, intind] = downsampling(nx, maxn, int_tol, M_kernel, ts_index, k_x);
    % Diffusion kernel - temperature, T
    [M_kernelR_T, indr_T, intind_T] = downsampling(nx, maxn, int_tol, M_kernel_T, ts_index_T, k_x);
    % Slip kernel - \delta_x
    [M_kernelR_deltax, indr_deltax, intind_deltax] = downsampling(nx, maxn, int_tol, M_kernel_deltax, ts_index_deltax, k_x);
    M_kernel_N = M_kernelR;
    M_kernel_S = M_kernelR;
    M_kernel_N_T = M_kernelR_T;
    M_kernel_S_T = M_kernelR_T;
    M_kernel_deltax_N = M_kernelR_deltax;
    M_kernel_deltax_S = M_kernelR_deltax;
    % Purge the memory for this matrix
    clearvars M_kernel
    clearvars M_kernel_deltax
    clearvars M_kernelR
    clearvars M_kernelR_deltax
    clearvars M_kernel_T
    clearvars M_kernelR_T
    maxM = zeros(1,length(k_x));
    for i_xx = 1:nx
    maxM(i_xx) = length(M_kernel_N{i_xx}(:, 1));
    end
    maxMT = zeros(1,length(k_x));
    for i_xx = 1:nx
    maxMT(i_xx) = length(M_kernel_N_T{i_xx}(:, 1));
    end
    maxMs = zeros(1,length(k_x));
    for i_xx = 1:nx
    maxMs(i_xx) = length(M_kernel_deltax_N{i_xx}(:, 1));
    end

    % Save the kernel to the .mat file
    kernel_list.ts_index = ts_index;
    kernel_list.ts_index_T = ts_index_T;
    kernel_list.ts_index_deltax = ts_index_deltax;
    kernel_list.indr = indr;
    kernel_list.indr_T = indr_T;
    kernel_list.indr_deltax = indr_deltax;
    kernel_list.M_kernel_N = M_kernel_N;
    kernel_list.M_kernel_S = M_kernel_S;
    kernel_list.M_kernel_N_T = M_kernel_N_T;
    kernel_list.M_kernel_S_T = M_kernel_S_T;
    kernel_list.M_kernel_deltax_S = M_kernel_deltax_S;
    kernel_list.M_kernel_deltax_N = M_kernel_deltax_N;
    kernel_list.maxM = maxM;
    kernel_list.maxMs = maxMs;
    kernel_list.maxMT = maxMT;  
    if nt_max < 10e4
        save(fullfile(kernel_folder, kernel_name), 'kernel_list');
        fprintf('Kernel saved to "%s".\n', fullfile(kernel_folder, kernel_name));
    else
        fprintf('Kernel size too large');
    end
else
    % If the folder exists, look for .mat files
    fprintf('Folder "%s" already exists.\n', kernel_folder);
    
    % Get the list of .mat files in the folder
    matFiles = dir(fullfile(kernel_folder, '*.mat'));
    
    if isempty(matFiles)
        fprintf('No .mat files found in the folder.\n');
    else
        fprintf('Found %d .mat file(s):\n', length(matFiles));
        for k = 1:length(matFiles)
            fprintf('  %s\n', matFiles(k).name);
            % Load the .mat file (example)
            load(fullfile(kernel_folder, matFiles(k).name));
            disp('Loaded data:');
            disp(kernel_name);
        end
        ts_index = kernel_list.ts_index;
        ts_index_T = kernel_list.ts_index_T;
        ts_index_deltax = kernel_list.ts_index_deltax;
        indr = kernel_list.indr;
        indr_T = kernel_list.indr_T;
        indr_deltax = kernel_list.indr_deltax;
        M_kernel_N = kernel_list.M_kernel_N;
        M_kernel_S = kernel_list.M_kernel_S;
        M_kernel_N_T = kernel_list.M_kernel_N_T;
        M_kernel_S_T = kernel_list.M_kernel_S_T;
        M_kernel_deltax_S = kernel_list.M_kernel_deltax_S;
        M_kernel_deltax_N = kernel_list.M_kernel_deltax_N;
        maxM = kernel_list.maxM;
        maxMs = kernel_list.maxMs;
        maxMT = kernel_list.maxMT;  
        clearvars kernel_list;
    end
end

jn_hat_storage_NM = zeros(nx, nt_max);
jn_hat_storage_SM = zeros(nx, nt_max); % no need to store everything
hat_delta_x_storage=zeros(nx, nt_max);
M_conv_S = zeros(size(k_x));
M_conv_N = zeros(size(k_x));
M_conv_delta_x = zeros(size(k_x));
% Convolution matrix for sigma_xz
M_conv_sxz_delta_x = zeros(size(k_x));
M_conv_sxz = zeros(size(k_x));
% Convolution matrix for sigma_zz
M_conv_szz = zeros(size(k_x));

jn_hat_storage_NM_T = zeros(nx, nt_max);
jn_hat_storage_SM_T = zeros(nx, nt_max);
M_conv_S_T = zeros(size(k_x));
M_conv_N_T = zeros(size(k_x));


ind_up = zeros(size(k_x)) + 1;
ind_up_s = zeros(size(k_x)) + 1;
ind_up_T = zeros(size(k_x)) + 1;
i_xv = 1:length(k_x);
i_xv_s = 1:length(k_x);
i_xv_T = 1:length(k_x);
disp("Convolution pre-allocation done")
%% ODE approximation for stress
load('coefficients.mat', 'a_opt', 's_i_opt');
% Define the reconstructed approximation function
reconstructed_approximation = @(s) ...
    sum(a_opt .* exp(-s ./ s_i_opt'), 1);
modelOrder = size(a_opt,2);
ro1 = a_opt;
damping1 = -1 ./ s_i_opt;
ODES2_sxz = zeros(nx, modelOrder); 
ODES2_sxz_delta_x = zeros(nx, modelOrder); 
ODES2_szz = zeros(nx, modelOrder);

disp("ODE approximation for convolution of stress done")

delta_T = 0;
i_tk = 1;
local_cumutime = 0;
% Knock on the sbi coupling
sbi_step = 1;
% Re-compute initial shear strain rate based on shear stress
gamma0 = 2 * gamma0 .* sinh(((sigma_xz + taub - G/(2*cs)*repmat(trapz(yy, gamma0, 2), 1, length(yy))) ./ (sigma_zz0 + sigma_zz - pa - p)) ./ (a - b)) ./...
                exp(fr0 ./ (a - b));
% Remove the boundary where the SBI boundary is resolved
gamma0 = gamma0(:,2:end-1);
% Enforce kinematic loading region

gamma          = zeros(nx, ny-2) + const_gamma0;
gamman         = zeros(nx, ny-2) + const_gamma0;
% Accumulated shear strain [1]
shear_strain = zeros(nx,ny-2);
yy_gammadot = yy(2:end-1);

% B. C.
% Periodic boundary for eastern and western sides
i1  = [1:nx]; % i
i2 = [nx 1:nx-1]; % i-1
i3 = [2:nx 1]; % i+1
% Grid indexing for y-direction
j = 2:ny-1;

dx2_inv = 1/dx^2;
dy2_inv = 1/dy^2;
rho_c_inv = 1 / rho_c;
alpha_th_dx2 = alpha_th * dx2_inv;
alpha_th_dy2 = alpha_th * dy2_inv;
dil_factor = eps_pl / beta;
% Temperature update
Te = @(i_t, Tn, gamman, tau) (tau.* gamman(i1,:) * rho_c_inv)+ ...
            alpha_th_dx2 .* ((Tn(i3,j)-Tn(i1,j)) - (Tn(i1,j)-Tn(i2,j)))+ ...
            alpha_th_dy2 .* ((Tn(i1,j+1)-Tn(i1,j)) - (Tn(i1,j)-Tn(i1,j-1)));
% Pore pressure update
pe = @(i_t, pn, gamman, gamma, delta_T,dtp) Lambda * delta_T  -  dil_factor .* (gamma(i1,:) - gamman(i1,:)) ./ (dtp.*gamman(i1,:))  + ...
            (M_alpha_hx(i3,j).*(pn(i3,j)-pn(i1,j)) - M_alpha_hx(i1,j).*(pn(i1,j)-pn(i2,j))) * dx2_inv + ...
            (M_alpha_hy(i1,j+1).*(pn(i1,j+1)-pn(i1,j)) - M_alpha_hy(i1,j).*(pn(i1,j)-pn(i1,j-1))) * dy2_inv;

% 1d
Te1d = @(i_t, Tn, gamman, tau, i_dsp) (tau.* gamman(i_dsp,:) / rho_c) + ...
            alpha_th .* ((Tn(i_dsp,j+1)-Tn(i_dsp,j)) - (Tn(i_dsp,j)-Tn(i_dsp,j-1))) / dy^2 ;

% Pore pressure update
pe1d = @(i_t, pn, gamman, gamma, delta_T,dtp, i_dsp) Lambda * delta_T  -  eps_pl * (gamma(i_dsp,:) - gamman(i_dsp,:)) ./ (dtp*beta.*gamman(i_dsp,:))  + ...
            (M_alpha_hy(i_dsp,j+1).*(pn(i_dsp,j+1)-pn(i_dsp,j)) - M_alpha_hy(i_dsp,j).*(pn(i_dsp,j)-pn(i_dsp,j-1))) / dy^2;

Intime_monitoring_rupture;

% Max iteration times for assuring linearization results converge
iter_max = 4;
iter_event = 0;
F_gammamax2_storage_plot = zeros(iter_max,1);
i_tF=0;
F_gammamax2 = zeros(1, ny - 2);
F_gamma =  zeros(nx, ny - 2);
penalty_gamma =  zeros(nx, ny - 2);
tic;
max_iter = 200; % Maximum number of iterations
tol = 1e-5;    % Convergence tolerance(taub - pa)
omega = 1;
i_t_unconverged = 0;
V = zeros(nx,1) + V0;
tauf = zeros(nx,1) + fr0*(sigma_zz0 - pa) - V0*eta;
maxit = 100;
sigma_xz_temp = zeros(1,nx);
sigma_zz_temp = zeros(1,nx);
sigma_xz_temp_nodeltax = zeros(1,nx);
tauhistg = sigma_xz_temp';
Vm = V';

nx_left = round(0.3*nx);
nx_right = round(0.7*nx);
E_bd = zeros(size(delta_x));

% Initialize time tracking
timev = [];
timev_post = [];

n_fig = 0;
inject_switch = 0;

dy_trap = diff(yy_gammadot(1:2)); % Assuming uniform spacing
inject_domain = [nx/2-2:nx/2+2];


%% Time stepping
%profile on
while cumutime < t_end   
    % FDM part
    i_t    = i_t + 1;

    % predictor 
    gammag = gamma;
    deltaxg = delta_x + dt * V';
    tauhistg = sigma_xz_temp_nodeltax + real(ifft2(akdhat4 .* fft2(deltaxg)));
    %pg = Inject_byrate(p, beta, t_end, i_t, dt, 5e-5, 1/2, 1/2, dx, dy, nx, ny, 3/nx, "Line");

    % previous time step
    pn     = p;
    Tn     = T;
    gamman = gamma;
    Vn = V;
    taufn = tauf;

    % Shear force balance without resolving T and p
    for i_x = 1:nx
        pint = p(i_x,j) + pa;
        [tauf_result, nb_it, delta_tauf] = ...
        solvestrength(tauf(i_x), taub0(i_x) + tauhistg(i_x), const_gamma0, A, sigma_zz0  + sigma_zz_temp(i_x), eta, pint, z, dy_trap, tol, maxit);  
        tauf(i_x) = tauf_result;
        if (i_x>=nx_left) && (i_x<=nx_right)
            V(i_x) = integratevelocity(z, const_gamma0, tauf(i_x), A, sigma_zz0  + sigma_zz_temp(i_x), pint, dy_trap);
        end
    end
    
    V(1:nx_left-1) = V0;
    V(nx_right+1:end) = V0;
    %
    deltaxg = delta_x + 0.5 * dt * (V + Vn)';
    tauhistg = sigma_xz_temp_nodeltax + real(ifft2(akdhat4 .* fft2(deltaxg)));
    gdot = zeros(1,ny-2);
    gamma = zeros(nx, ny-2);
    
    for i_x = 1:nx
        pint = p(i_x,j) + pa;
        [tauf_result, nb_it, delta_tauf] = ...
        solvestrength(tauf(i_x), taub0(i_x) + tauhistg(i_x), const_gamma0, A, sigma_zz0  + sigma_zz_temp(i_x), eta, pint, z, dy_trap, tol, maxit); 
        tauf(i_x) = tauf_result;
        if (i_x>=nx_left) && (i_x<=nx_right)
            %V(i_x) = integratevelocity(z, const_gamma0, tauf(i_x), A, sigma_zz0 + sigma_zz_temp(i_x), pint, dy_trap);
            gdot = strainrate(z, const_gamma0, tauf(i_x), A, sigma_zz0  + sigma_zz_temp(i_x), pint);
            %gdot = gdot ./ fast_trapz(gdot, dy_trap) .* V(i_x);
            V(i_x) = fast_trapz(gdot, dy_trap);
            gamma(i_x,:) = gdot;
        end
        %nb_it_(i_x) = nb_it;
    end
    
    %size(tauf_)
    gamma(1:nx_left-1,:) = const_gamma0; gamma(nx_right+1:end,:) = const_gamma0;
    V(1:nx_left-1) = V0; V(nx_right+1:end) = V0;
    tauf_ = taub(:) + tauhistg(:) - eta * V(:);
    FV = Vgamma(V, z, repmat(tauf,1,ny-2) , sigma_zz0 + sigma_zz(:,j) - p(:,j) - pa, parser.a, parser.b, V0/h, dy_trap);
    FV = FV(nx_left:nx_right);

    Vm = 0.5 * (Vn + V)';
    delta_x    = delta_x + Vm*dt;
    delta_x_incre = Vm*dt;
    
    % predictor
    T(i1, j) = Tn(i1,j) + Te(i_t*dt, Tn, tauf, gamma) * dt;
    delta_T  = (T(i1,j) - Tn(i1,j))/dt;
    p(i1, j) = pn(i1, j) + pe(i_t*dt, pn, gamman, gamma, delta_T, dt)*dt;
    
    % T and p at midpoint
    p_mid    = (p + pn)*0.5;
    T_mid    = (T + Tn)*0.5;
    gamma_mid = (gamma + gamman)*0.5;
    
    % T and p at the next half time step 
    T(i1, j) = T_mid(i1,j) + Te(i_t*dt, Tn, tauf, gamma) * 0.5 * dt;
    delta_T  = (T(i1,j) - T_mid(i1,j))/(0.5*dt);
    p(i1, j) = p_mid(i1, j) + pe(i_t*dt, pn, gamma_mid, gamma, delta_T, 0.5*dt)* 0.5 *dt;
    
    % Injection line source
    %p = Inject_byrate(p, beta, 0.0050, i_t, dt, 5e-5, 1/2, 1/2, dx, dy, nx, ny, 3/nx, "Line");

    %p(inject_domain,2:end-1) = Inject_bypressure(cumutime+dt, t_end*0.75, t_end, 1, 110, 20);
    shear_strain = shear_strain + gamma*dt;
    
    % SBI part
    % profile on
    if sbi_step 
    
    % Cal jn in the coupled interface
    jn_N = -M_kappa_y(2:end,end-1) .* (p(:,end) - p(:,end-1))/dy;
    jn_S = -M_kappa_y(2:end,2).* (p(:,2) - p(:,1))/dy;
    if i_tk>1
    jn_hat_preN = jn_hat_storage_NM(:, i_tk-1);
    jn_hat_preS = jn_hat_storage_SM(:, i_tk-1);
    
    else
    jn_hat_preN = 0;
    jn_hat_preS = 0;   
    end
    jn_hat_N = fft2(jn_N);
    jn_hat_S = fft2(jn_S);
    jn_hat_storage_NM(:, i_tk) = 0.5*(jn_hat_N + jn_hat_preN);  % integration scheme change
    jn_hat_storage_SM(:, i_tk) = 0.5*(jn_hat_S + jn_hat_preS);  
    hat_delta_x = fft2(delta_x);
    hat_delta_x_storage(:, i_tk) = hat_delta_x;
    % Convolution for computing pressure
    Iup = (i_tk - ind_up)./ts_index > tol_up;       
    i_xred = i_xv(Iup);                           
    for i_x = i_xred                              
            % if i_t == 1 || (i_t - ind_up(i_x))/ts_index(i_x) > tol_up
            if i_tk < ts_index(1,i_x)
                indvT = flip(indr{i_x}(:,1));
                Inonneg = -indvT + i_tk > 0;
                indvT = indvT(Inonneg);
                M_conv_N(i_x) = jn_hat_storage_NM(i_x, indvT) * ((M_kernel_N{i_x}((maxM(i_x) - length(indvT)+1):maxM(i_x))));
                M_conv_S(i_x) = jn_hat_storage_SM(i_x, indvT) * ((M_kernel_S{i_x}((maxM(i_x) - length(indvT)+1):maxM(i_x))));
                M_conv_delta_x(i_x) = hat_delta_x_storage(i_x, indvT) * ((M_kernel_N{i_x}((maxM(i_x) - length(indvT)+1):maxM(i_x))));
            else
                indvT = i_tk  -  indr{i_x}(:,1) + 1 ; % pick up the nodes    
                M_conv_N(i_x) = jn_hat_storage_NM(i_x,indvT) * ( M_kernel_N{i_x}  );
                M_conv_S(i_x) = jn_hat_storage_SM(i_x,indvT) * ( M_kernel_S{i_x}  );
                M_conv_delta_x(i_x) = hat_delta_x_storage(i_x,indvT) * ( M_kernel_N{i_x}  );
            end
            ind_up(i_x) = i_tk;
    end
    % p - SBI boundary 
    % Northern 
    %delta_x = zeros(nx,1)'; delta_x(nx/4:3*nx/4) = normpdf(xx(nx/4:1:3*nx/4),0,30); hat_delta_x = fft2(delta_x);
    p(:,end) = real(ifft2(-akdhat1.* hat_delta_x  ...
               -akdhat2 .* dt0*kernelfac.*M_conv_delta_x  ...
               -akdhat3 .* dt0*kernelfac.*M_conv_N)); 
    % Separting the contribution of poroelasticity in p
    p_poro_end = p(:,end)' - real(ifft2(akdhat3 .* dt0*kernelfac*M_conv_N));
    p_nodeltax_end = real(ifft2(-akdhat2 .* dt0*kernelfac.*M_conv_delta_x  ...
               -akdhat3 .* dt0*kernelfac.*M_conv_N));
    %p_akdhat1N = real(ifft2(-akdhat1.* hat_delta_x));
    %p_akdhat2N = real(ifft2(-akdhat2 .* dt0*kernelfac.*M_conv_delta_x));
    %p_akdhat3N = real(ifft2(-akdhat3 .* dt0*kernelfac.*M_conv_N));
    % Southern
    p(:,1)   = real(ifft2(akdhat1.* hat_delta_x  ...
               +akdhat2 .* dt0*kernelfac.*M_conv_delta_x  ...
               +akdhat3 .* dt0*kernelfac.*M_conv_S));
    %p_akdhat1S = real(ifft2(+akdhat1.* hat_delta_x));
    %p_akdhat2S = real(ifft2(+akdhat2 .* dt0*kernelfac.*M_conv_delta_x));
    %p_akdhat3S = real(ifft2(+akdhat3 .* dt0*kernelfac.*M_conv_S));
    % Separting the contribution of poroelasticity in p
    p_poro_first = p(:,1)' - real(ifft2(akdhat3 .* dt0*kernelfac*M_conv_S));
    p_nodeltax_first = real(ifft2(akdhat2 .* dt0*kernelfac.*M_conv_delta_x  ...
               +akdhat3 .* dt0*kernelfac.*M_conv_S));

    % ODE approximation
    ODES2_sxz = (ODES2_sxz + ro1 .* ((jn_hat_N + jn_hat_S) .* dt0)) ./ ...
                 (1 - alpha_hostrock .* damping1 .* k_x'.^2 .* dt0);
    ODES2_sxz_delta_x = (ODES2_sxz_delta_x + ro1 .* hat_delta_x' .* dt0) ./ ...
                         (1 - alpha_hostrock.*damping1.*k_x'.^2.*dt0);
    ODES2_szz = (ODES2_szz + ro1 .* (jn_hat_N - jn_hat_S) .* dt0) ./ ...
                 (1 - alpha_hostrock.*damping1.*k_x'.^2.*dt0);
    M_conv_sxz = sum(ODES2_sxz'); 
    M_conv_sxz_delta_x = sum(ODES2_sxz_delta_x');
    M_conv_szz = sum(ODES2_szz');

    % sigma_xz - SBI boundary 
    sigma_xz_temp = real(ifft2(akdhat4.* hat_delta_x +  akdhat5 .* kernelfac .* M_conv_sxz_delta_x ...
                      + akdhat6 .* kernelfac .* M_conv_sxz));
    sigma_xz_temp5 = real(ifft2(akdhat5 .* kernelfac .* M_conv_sxz_delta_x));
    sigma_xz_temp6 = real(ifft2(akdhat6 .* kernelfac .* M_conv_sxz));
    sigma_xz_temp_nodeltax = sigma_xz_temp - real(ifft2(akdhat4.* hat_delta_x));
    % constant sigma_xz across the fault
    sigma_xz = repmat(sigma_xz_temp', 1, ny);
    % sigma_zz - FD boundary - Periodic by SBIM solutions
    sigma_zz_temp = real(ifft2(akdhat7 .* kernelfac .* M_conv_szz));
    % constant sigma_zz across the fault
    sigma_zz = repmat(sigma_zz_temp', 1, ny);

    % T - zero gradient Neumann, in Platt, they assum the heat can flux into the bulk
    jn_N_T = -M_kappa_T_y(2:end,end-1) .* (T(:,end) - T(:,end-1))/dy;
    jn_S_T = -M_kappa_T_y(2:end,2).* (T(:,2) -T(:,1))/dy;
    if i_tk>1
    jn_hat_preN_T = jn_hat_storage_NM_T(:, i_tk-1);
    jn_hat_preS_T = jn_hat_storage_SM_T(:, i_tk-1); 
    %diff_dt = (timev(end) - i_tk*dt0); 
    else
    jn_hat_preN_T = 0;
    jn_hat_preS_T = 0;   
    end
    jn_hat_storage_NM_T(:, i_tk) = 0.5*(fft2(jn_N_T) + jn_hat_preN_T);  
    jn_hat_storage_SM_T(:, i_tk) = 0.5*(fft2(jn_S_T) + jn_hat_preS_T);  
    
    % Convolution for computing pressure
    Iup_T = (i_tk - ind_up_T)./ts_index_T > tol_up;       
    i_xred_T = i_xv_T(Iup_T);                           
    for i_x = i_xred_T                             
            if i_tk < ts_index_T
                indvT_T = flip(indr_T{i_x}(:,1));
                Inonneg = -indvT_T + i_tk > 0;
                indvT_T = indvT_T(Inonneg);
                M_conv_N_T(i_x) = jn_hat_storage_NM_T(i_x, indvT_T) * ((M_kernel_N_T{i_x}((maxMT(i_x) - length(indvT_T)+1):maxMT(i_x))));
                M_conv_S_T(i_x) = jn_hat_storage_SM_T(i_x, indvT_T) * ((M_kernel_S_T{i_x}((maxMT(i_x) - length(indvT_T)+1):maxMT(i_x))));
            else
                indvT_T = i_tk  -  indr_T{i_x}(:,1) + 1 ; % pick up the nodes    
                M_conv_N_T(i_x) = jn_hat_storage_NM_T(i_x,indvT_T) * M_kernel_N_T{i_x};
                M_conv_S_T(i_x) = jn_hat_storage_SM_T(i_x,indvT_T) * M_kernel_S_T{i_x};
            end
            ind_up_T(i_x) = i_tk;
    end    

    T(:,end) = real(ifft2(-akdhatT*dt0*M_conv_N_T));    
    T(:,1) = real(ifft2(akdhatT*dt0*M_conv_S_T));
    i_tk = i_tk + 1;

    local_cumutime = 0;

    else
    
    % undrained poroelastic response
    p(:,end) = p_nodeltax_end + real(ifft2(-akdhat1.* hat_delta_x)); 
    p(:,1)   = p_nodeltax_first + real(ifft2(akdhat1.* hat_delta_x));
    sigma_xz_temp = sigma_xz_temp_nodeltax + real(ifft2(akdhat4.* hat_delta_x));
   
    % constant sigma_xz across the fault
    sigma_xz = repmat(sigma_xz_temp', 1, ny);
    end
    % profile off
    % profile viewer
    
    dtp = dt;
    dt = dt0;

    if (max(FV) > 1 / 128) && (dt / 2 > dt0 / 256)
        dt = 0.5 * dt;
    else
        dt = dt0;
    end

    if (local_cumutime + dt) == dt0
        sbi_step = 1;
    elseif (local_cumutime + dt) < dt0
        sbi_step = 0;
    else
        dt = dt0 - local_cumutime;
        sbi_step = 1;
    end

    local_cumutime = local_cumutime + dt;

    if  any(isnan(gamma(:)) | isinf(gamma(:)))
        disp("Model on the verge of explosion");
        disp("End simulations");
        break;
    else
        shear_strain_final_plot = shear_strain + (gamma + gamman)*dt;
    end
    %% Results dodged for post-processing    
    cumutime = cumutime + dt;

    % time
    timev(i_t) = cumutime;
    % Save variables recurrently for post-processing
    
    % breakdown energy
    E_bd = E_bd + ((taufn + tauf) .* 0.5 - taub)' .* delta_x_incre;
    % loglog(delta_x_storage_plot(nx/2,:), abs(E_bd_plot(nx/2,:)))
    if rem(i_t,400)==0
        
        % count fig
        n_fig = n_fig + 1;

        % t saved for postprocessing
        writeout.timev_post(n_fig) = cumutime;

        % Slip velocity profile [m/s]
        writeout. V_storage_plot(:,n_fig) = V;

        % Slip distance profile [m]
        writeout.delta_x_storage_plot(:,n_fig) = delta_x;
        
        % Max gamma
        writeout.gamma_max_plot(:, n_fig) = max(gamma, [], 2);

        % Shear strain rate - section positioned at x = 3/4*L
        writeout.gamma_storage_plot_2o3(:, n_fig) = gamma(round(nx*0.6),:);

        % Shear strain rate - section positioned at x = 1/2*L
        writeout.gamma_storage_plot_1o2(:, n_fig) = gamma(nx/2,:);

        % Shear strain rate - section positioned at x = 1/4*L
        writeout.gamma_storage_plot_1o3(:, n_fig) = gamma(round(nx*0.4),:);

        % Shear stress on the fault
        tau_shear = tauf;
        writeout.tau_shear_storage(:, n_fig) = tauf;
        sigma_xz5_storage(:, n_fig) = sigma_xz_temp5;
        sigma_xz6_storage(:, n_fig) = sigma_xz_temp6;
        % Normal stress on the fault
        sigma_normal = sigma_zz0+sigma_zz(:,2:end-1)-pa-p(:,2:end-1);
        writeout.sigma_normal_storage(:, n_fig) = mean(sigma_normal,2);

        % Friction check - if negative friction is generated
        writeout.fric_storage_plot(:,n_fig)=min((tauf ./ sigma_normal)');
        
        % Pore pressure - down boundary
        writeout.p_south_storage_plot(:, n_fig) = p(:,1);

        % Pore pressure contributed by poroelasticity - down boundary
        writeout.p_poro_south_storage_plot(:, n_fig) = p_poro_first;

        % Pore pressure - up boundary
        writeout.p_north_storage_plot(:, n_fig) = p(:,end);

        % Pore pressure contributed by poroelasticity - up boundary
        writeout.p_poro_north_storage_plot(:, n_fig) = p_poro_end;

        % Pore pressure - section positioned at x = 1/2*L
        writeout.p_storage_plot(:, n_fig) = p(nx/2,:);

        % Temperature - section positioned at x = 1/2*L
        writeout.T_storage_plot(:, n_fig) = T(nx/2,:);

        % Shear stress contribution computed by SBIM
        writeout.sigma_xz_storage_plot(:, n_fig) = sigma_xz_temp;
        
        % Normal stress contribution computed by SBIM
        writeout.sigma_zz_storage_plot(:, n_fig) = sigma_zz_temp;
        
        % Characteristic pore pressure on the fault
        writeout.p_character(:,n_fig) = sigma_zz0(:) + sigma_zz_temp(:) - tauf(:)./(A*log(V(:)./V0) + fr0)-pa;
        
        % Log weighted
        writeout.p_weight_loggamma(:,n_fig) = sum(p(:,2:end-1).*log(gamma./gamma_c + exp(-5)),2)./sum(log(gamma./gamma_c + exp(-5)),2);
        
        % Log weighted
        asinhgamma = asinh(gamma./(2*gamma_c).*exp(fr0./A));
        writeout.p_weight_loggamma_regularized(:,n_fig) = sum(p(:,2:end-1).* asinhgamma,2) ./sum(asinhgamma,2);

        % Mean pore pressure on the fault
        writeout.p_mean(:,n_fig) = mean(p,2);
        
        % Shear strain rate weighted pore pressure on the fault
        writeout.p_weight_gamma(:,n_fig) = sum(p(:,2:end-1).*gamma,2)./sum(gamma,2);
        
        % Max pore pressure on the fault
        writeout.p_max(:,n_fig) = max(p');
        
        % breakdown energys
        writeout.E_bd_plot(:, n_fig) = E_bd;

        % strain 
        writeout.max_strain_plot(:, n_fig) = max(shear_strain);
        % Terminate simulation on the edge of explosion
        
        % Accumulated shear strain
        if rem(i_t,freq_post)==0
            
            plot_physics_field; 
            drawnow;
            frame = getframe(fig); % Capture the frame
            writeVideo(videoFile2, frame); % Write the frame to the video

            % dump fields
            p_dump = sprintf('%s_%d.mat', 'p', i_t/freq_post);
            T_dump = sprintf('%s_%d.mat', 'T', i_t/freq_post);
            f_dump = sprintf('%s_%d.mat', 'f', i_t/freq_post);
            
            save(fullfile(parentFolder, 'post', outputFolder, p_dump), 'p');
            save(fullfile(parentFolder, 'post', outputFolder, T_dump), 'T');
            save(fullfile(parentFolder, 'post', outputFolder, f_dump), 'friction');
            if rem(i_t,freq_post*50) == 0
                save(fullfile(parentFolder, 'post', outputFolder, 'AnisoPermea.mat'), 'writeout');
            end
        end
        % Visualization
        
    end
 
 
end
%profile off
%profile viewer
% Consumed time
writeout.t_wallclock = toc;
save(fullfile(parentFolder, 'post', outputFolder, 'AnisoPermea.mat'), 'writeout');
%% Save output

% Save the table as a mat file

close(videoFile2);
%close all; 
disp('End Simulations');
