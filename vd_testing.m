%% Velocity Direction Testing
% examing bug with velocity direction/angular momentum vector direction
clc 
clear
close all

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

%% Lambert set up and solve
    % Create Random Initial Conditions


    [r1, v1, r2, v2, delta_t] = lambert_conditions(mu);
    [a_v, e_v, p_v, i_v, Omega_1_v, ~, ~] = get_oe(r1, v1, mu);
    [i_1, Omega_11, Omega_21] = newton_angles(r1, r2, a_v, e_v, J_2, mu, alpha, delta_t);
    
    T = 2*pi*sqrt((a_v^3)/mu);
    ratio = delta_t/T;
    
    % Lambert Solve
    N = 15;
    [a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);
    error = norm(v1_L - v1)/norm(v1);


 %% Angular Momentum Plot
 h1_v = cross(r1,v1)/norm(cross(r1,v1));
 h2_v = cross(r2,v2)/norm(cross(r2,v2));

 h1_L =[sin(i_1)*sin(Omega_11); -sin(i_1)*cos(Omega_11); cos(i_1)];
 test = cross(r1,r2)/norm(cross(r1,r2));
    
 % Plot
 figure
 plot3([0 h1_v(1)], [0 h1_v(2)], [0, h1_v(3)], 'r-')
 hold on
 plot3([0 h2_v(1)], [0 h2_v(2)], [0, h2_v(3)], 'y-')
 plot3([0 test(1)], [0 test(2)], [0, test(3)], 'g--')
 plot3([0 h1_L(1)], [0 h1_L(2)], [0, h1_L(3)], 'b--')
