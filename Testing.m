%% This Scipt is for testing the Generalized Labert's Problem 
clc;
clear
close all
%% Given

% time of flight
delta_t = 3000;        % [s]   

% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

% Position vectors
r1 = [7000; 50; 300];    % [km] 

%% Test Case
% Initial Velocity
v1 = [1; -8; -2];

% Kepler Solve
[r2,v2] = pkepler(r1, v1, delta_t, 0, 0);


theta = acos(dot(r1/norm(r1), r2/norm(r2)));

% Calculate Velocity Magnitudes
v1_mag_real = norm(v1);
v2_mag_real = norm(v2);

% Inital and Final RAAN and inclination
[a_v, e_v, p_v, i_v, Omega_1_v, ~, ~] = get_oe(r1, v1, mu);
[~, ~, ~, ~, Omega_2_v, ~, ~] = get_oe(r2, v2, mu);
[p,a_p,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r1,v1,mu);
T = 2*pi*sqrt(a_v^3/mu);
ratio = delta_t/T;

[a, v1_L, v2_L] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, 100);
[a_L, ~, ~, i_v, Omega_1_l, ~, ~] = get_oe(r1, v1_L', mu);
[i_1, Omega_11, Omega_21] = newton_angles(r1, r2, a_v, e_v, J_2, mu, alpha, delta_t);
[i_2, Omega_12, Omega_22] = newton_angles(r1, r2, a, e, J_2, mu, alpha, delta_t);
%% Test things
% cord length
c = norm(r1 - r2);
% semi-perimeter
semi = .5*(norm(r1) + norm(r2) + c);
[psi, cphi] = lagrange_angles(semi,c,a_v);
e_test = get_eccentricity(r1, r2, theta, psi, a_v);
p_test = norm(r1)*norm(r2)*sin(theta/2)^2/(a*sin(psi)^2);

e = get_eccentricity(r1, r2, theta, psi, a_v);
[i, Omega_1_L, Omega_2_L] = newton_angles(r1, r2, a, e, J_2, mu, alpha, delta_t);

%% Results
% Error in semi major axis
a_error = abs(a_v - a)/a_v;
% Error in velocity magnitude
v1_error = abs(v1_mag_real - v1_L);
v2_error = abs(v2_mag_real - v2_L);





