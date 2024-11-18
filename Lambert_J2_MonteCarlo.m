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

% Vector Direction
v1_direction = rand(1,3);
v1_direction = v1_direction/norm(v1_direction);
r1_direction = rand(1,3);
r1_direction = r1_direction/norm(r1_direction);

% Vector Magnitude
% Maximum radius set to GEO radius
r_max = 42164;
r1_mag = alpha + (r_max - alpha)*rand;

% Maximum velocity set to leo escape velocity
e_max = .5;
v_max = sqrt(mu*(1 + e_max)/r1_mag);
v_circle = sqrt(mu/r1_mag);
v1_mag = v_circle + (v_max - v_circle)*rand();

% Final Vectors
r1 = r1_mag*r1_direction;
v1 = v1_mag*v1_direction;

% Random Time of Flight
tof_min = 500;
tof_max = 2*pi*sqrt(r_max^3/mu);
delta_t = tof_min + (tof_max - tof_min)*rand();

%% Find Final Condtions
[r2,v2] = pkepler(r1, v1, delta_t, 0, 0);
v2_mag = norm(v2);

[a_v, e_v, p_v, i_v, Omega_1_v, ~, ~] = get_oe(r1, v1, mu);
%% Lambert Solve
N = 10;
[a_L, v1_L, v2_L] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, N);

%% Errors
error_a_abs = abs(a_L - a_v);
error_v1_abs = abs(v1_L - v1_mag);
error_v2_abs = abs(v2_L - v2_mag);