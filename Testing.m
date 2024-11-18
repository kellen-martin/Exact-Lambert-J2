%% This Scipt is for testing the Generalized Labert's Problem 
clc;
clear
close all
%% Given

% time of flight
delta_t = 500;        % [s]   

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
v1 = [-1; 7.5;-5];

% Kepler Solve
[r2,v2] = pkepler(r1, v1, delta_t, 0, 0);


theta = acos(dot(r1/norm(r1), r2/norm(r2)));

% Calculate Velocity Magnitudes
v1_mag_real = norm(v1);
v2_mag_real = norm(v2);

% Inital and Final RAAN and inclination
[a_v, e_v, p_v, i_v, Omega_1_v, ~, ~] = get_oe(r1, v1, mu);
[~, ~, ~, ~, Omega_2_v, ~, ~] = get_oe(r2, v2, mu);
T = 2*pi*sqrt(a_v^3/mu);


[a, v1, v2] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, 100);


%% Results

% Absolute Error in velocity magnitude
v1_error = abs(v1_mag_real - v1);
v2_error = abs(v2_mag_real - v2);





