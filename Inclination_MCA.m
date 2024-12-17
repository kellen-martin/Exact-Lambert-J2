%% Monte Carlo Analysis
% Compares error and computation time of using the inclination at a0, (the
% center of the circular contour) with calculating inclination at each
% evaluation of f
clc
clear

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

%% Monte Carlo
% number of quaderature points
N = 64;
% number of trials
trials = 1000;


error_v1_rel_1newt = zeros(1,trials);
error_v2_rel_1newt = zeros(1,trials);
error_a_rel_1newt = zeros(1,trials);
solve_time_1newt = zeros(1,trials);

error_v1_rel = zeros(1,trials);
error_v2_rel = zeros(1,trials);
error_a_rel = zeros(1,trials);
solve_time = zeros(1,trials);

    for i=1:trials
    % Create Random Initial Conditions
    
    [r1, v1, r2, v2, delta_t] = lambert_conditions(mu);
    [a_v, ~, ~, ~, ~, ~, ~] = get_oe(r1, v1, mu);

    % Lambert Solve
    tic
    [a_1newt, v1_L_1newt, v2_L_1newt] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);
    solve_time_1newt(i) = toc;

    tic
    [a, v1_L, v2_L] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, N);
    solve_time(i) = toc;
    
    % Errors
    error_v1_rel_1newt(i) = norm(v1_L_1newt - v1)/norm(v1);
    error_v2_rel_1newt(i) = norm(v2_L_1newt - v2)/norm(v2);
    error_a_rel_1newt(i) = abs(a_1newt - a_v)/abs(a_v);

    error_v1_rel(i) = norm(v1_L - v1)/norm(v1);
    error_v2_rel(i) = norm(v2_L - v2)/norm(v2);
    error_a_rel(i) = abs(a - a_v)/abs(a_v);
    end

% Mean Error
a_error_mean_1newt = mean(error_a_rel_1newt);
v1_error_mean_1newt = mean(error_v1_rel_1newt);
v2_error_mean_1newt = mean(error_v2_rel_1newt);
solve_time_mean_1newt = mean(solve_time_1newt);
std_a_error_1newt = std(error_a_rel_1newt);

a_error_mean = mean(error_a_rel);
v1_error_mean = mean(error_v1_rel);
v2_error_mean = mean(error_v2_rel);
solve_time_mean = mean(solve_time);
std_a_error = std(error_a_rel);

%% Compare Results
a_error_diff = a_error_mean - a_error_mean_1newt;
v1_error_diff = v1_error_mean - v1_error_mean_1newt;
v2_error_diff = v2_error_mean - v2_error_mean_1newt;

comp_time_diff = solve_time_mean - solve_time_mean_1newt;
comp_time_reduction = comp_time_diff/solve_time_mean;
