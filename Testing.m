%% This Scipt is for testing the Generalized Labert's Problem 
clc;
clear
close all

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

% Time step
t_step = 1;        % [s]

%% Generate random orbit 
N = 64;
% [r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);
r1 = [5.759301829612105e+03;-5.824374577812693e+03;1.959376242451960e+03];
r2 = [7.267147194508008e+03;1.016451790376569e+04;-3.552235887329285e+03];
v1 = [6.663921814079755;2.872583327035716;-1.036325223281615];
v2 = [-3.945046393547707;2.099552425014366;-0.689536423412820];
delta_t = 3.959601521502576e+03;
[a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);


[r2_sim,v2_sim] = pkepler(r1, v1_L, delta_t, 0, 0);

error = norm(r2 - r2_sim)/norm(r2);


%% 



