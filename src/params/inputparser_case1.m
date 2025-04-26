% This parameter set as a reference 
% Geometry and grid
parser.L = 500; % Fault length and x-axis grid counts
parser.nx = 2^12; % Counts of the x-axis node
parser.ny = 150; % Counts of the y-axis node
parser.h = 0.0005; % Model thickness

% Poroelastic boundary 
% Spectral boundary integral
% Source: Wang, 2001; Noda, 2022
parser.vu = 0.34; % Undrained Poisson's ratio
parser.v = 0.25; % Drained Poisson's ratio
parser.G = 30e9; % Shear modulus
parser.cs = 3000; % Shear wave velocity
parser.alphaBiot = 0.47;
parser.B = ((1 + parser.vu) / (1 - 2 * parser.vu) - (1 + parser.v) / (1 - 2 * parser.v)) ...
           / ( parser.alphaBiot * (1 + parser.vu) / (1 - 2 * parser.vu)); % SkemptonBs coefficient
           
% 1 km Platt 2014
%parser.Lambda = 0.068e6;
%parser.rho_c = 2.7e6;
%parser.alpha_th = 0.7e-6;
%parser.alpha_hx = 7.15e-6;
%parser.alpha_hy = 7.15e-6;
%parser.eps_pl = 1.7e-4*0; % 
%parser.beta = 4.39e-10; % Compressibility (1/Pa) 

% 7 km Platt 2014
parser.Lambda = 0.30e6;
parser.rho_c = 2.7e6;
parser.alpha_th = 0.54e-6;
parser.alpha_hx = 6.71e-6;
parser.alpha_hy = 6.71e-6;
parser.eps_pl = 1.7e-4; % 
parser.beta = 2.97e-10; % Compressibility (1/Pa) 
parser.fr0 = 0.6;
parser.fr = 0.6;
parser.beta_T = 1e-6; 
parser.a = (0.040); % Rate-dependent
parser.b = (0.040 - 0.025); % Rate-dependent

parser.alpha_hostrock = parser.alpha_hy;

% Precomputed coefficients for implementation
A = parser.a - parser.b;
z = parser.fr0 / A;

% Time stepping
parser.endingtime = 0.075;
parser.sigma = 1/2^4;