%% Monte Carlo Analysis
clc
clear

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

%% Create Random Initial Conditions

[r1, v1, r2, v2, delta_t] = lambert_conditions(mu);
[a_v, e_v, p_v, i_v, Omega_1_v, ~, ~] = get_oe(r1, v1, mu);
[i_1, Omega_11, Omega_21] = newton_angles(r1, r2, a_v, e_v, J_2, mu, alpha, delta_t);

T = 2*pi*sqrt((a_v^3)/mu);
ratio = delta_t/T;

%% Test Things
theta = 2*pi - acos(dot(r1/norm(r1), r2/norm(r2)));
% cord length
c = norm(r1 - r2);
% semi-perimeter
semi = .5*(norm(r1) + norm(r2) + c);
[psi, cphi] = lagrange_angles(semi,c,a_v);
e_test = get_eccentricity(r1, r2, theta, psi, a_v);
p_test = norm(r1)*norm(r2)*sin(theta/2)^2/(a_v*sin(psi)^2);

%% Lambert Solve
N = 100;
[a, v1_mag, v2_mag] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, N);

%% Errors
v1_mag_real = norm(v1);
v2_mag_real = norm(v2);
error_v1_mag = abs(v1_mag - v1_mag_real)/v1_mag_real;
error_v2_mag = abs(v2_mag - v2_mag_real)/v2_mag_real;
error_a_rel = abs(a - a_v)/a_v;
