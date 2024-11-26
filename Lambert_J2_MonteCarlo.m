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

%% Lambert Solve
N = 100;
[a_o, ~, ~] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, N);

%% Errors
%error_a_abs = abs(a_L - a_v);
%error_v1_abs = abs(v1_L - v1_mag);
%error_v2_abs = abs(v2_L - v2_mag);