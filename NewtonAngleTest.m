%% Test of Newton Angles
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
check = true;
n = 0;
N = 128;


[r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);

% Orbital Elements at 1 and 2
[a1, e1, p1, i1, Omega1, omega1, f1] = get_oe(r1, v1, mu);
[a2, e2, p2, i2, Omega2, omega2, f2] = get_oe(r2, v2, mu);

% True Angular Momentum Vectors
h1 = cross(r1, v1);
h2 = cross(r2, v2);
h1_hat = h1/norm(h1);
h2_hat = h2/norm(h2);
r_cross_r = cross(r1, r2)/norm(cross(r1, r2));

% Find Orbital Planes
[in, Omega_1, Omega_2] = newton_angles(r1, r2, a1, e1, J_2, mu, alpha, delta_t);
hn1 = [sin(in)*sin(Omega_1), -sin(in)*cos(Omega_1), cos(in)];
% hn1_hat = hn1/norm(hn1);
hn2 = [sin(in)*sin(Omega_2), -sin(in)*cos(Omega_2), cos(in)];
% hn2_hat = hn2/norm(hn2);

val1 = dot(hn1, r1);
val_11 = hn1(1)*r1(1) + hn1(2)*r1(2) + hn1(3)*r1(3);
val2 = dot(hn2, r2');

x = [in, Omega_1, Omega_2];
    Beta = 3/2*((J_2*sqrt(mu)*alpha^2)/(a1^(7/2)*(1 - e1^2)^2));
F = @(x) [sin(x(1))*sin(x(2))*r1(1) - sin(x(1))*cos(x(2))*r1(2) + cos(x(1))*r1(3);
              sin(x(1))*sin(x(3))*r2(1) - sin(x(1))*cos(x(3))*r2(2) + cos(x(1))*r2(3);
              x(3) + Beta*cos(x(1))*delta_t - x(2)];
test = F(x);
 % Plot
 figure
 plot3([0 h1_hat(1)], [0 h1_hat(2)], [0, h1_hat(3)], 'r-')
 hold on
 plot3([0 hn1(1)], [0 hn1(2)], [0, hn1(3)], 'b-')

% [in, Omega_1, Omega_2] = newton_angles(r1, r2, a1, e1, J_2, mu, alpha, delta_t);


% Case 1 component 3 is flipped
% i = i - pi

% Case 2 component 1 & 2 are flipped
% i = -i