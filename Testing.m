%% This Scipt is for testing the Generalized Labert's Problem 
clc;
clear
close all
%% Given

% time of flight
delta_t = 2000;        % [s]   

% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

% Position vectors
r1 = [7000; 500; 30];    % [km] 

%% Test Case
% Initial Velocity
v1 = [-1; 7.5; 2];

% Kepler Solve
[r2,v2] = pkepler(r1, v1, delta_t, 0, 0);


% Calculate Velocity Magnitudes
v1_mag_real = norm(v1);
v2_mag_real = norm(v2);

% Inital and Final RAAN and inclination
[a_v, e_v, p_v, i_v, Omega_1_v, ~, ~] = get_oe(r1, v1, mu);
[~, ~, ~, ~, Omega_2_v, ~, ~] = get_oe(r2, v2, mu);



[a, v1, v2] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, 10);



%% Results

% Absolute Error in velocity magnitude
v1_error = abs(v1_mag_real - v1);
v2_error = abs(v2_mag_real - v2);





