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



%% Compute other quantites
% angle between r1 & r2
theta = acos(dot(r1/norm(r1), r2/norm(r2)));

% cord length
c = norm(r1 - r2);

% semi perimeter
semi = .5*(norm(r1) + norm(r2) + c);

%% More Tests
[i_test, Omega_1_test, Omega_2_test] = newton_angles_2(r1, r2, a_v, e_v, J_2, mu, alpha, delta_t);
[psi_test, c_phi_test] = get_angles(semi,c,a_v);

p_test = norm(r1)*norm(r2)*sin(theta/2)^2/(a_v*sin(psi_test)^2);
e_test = get_eccentricity(r1, r2, theta, psi_test, a_v);

%% Compute contour
a_min = c/2;
e_max = .8;
r_p = min(norm(r1), norm(r2));
r_a = -r_p*(e_max + 1)/(e_max - 1);
a_max = (r_a + r_p)/2;
Cr = (a_max - a_min)/2;
a0 = a_min + Cr;

%% Trapaziod Rule

N = 10; % [number of points]
sum1 = 0;
sum2 = 0;
for j=1:(N-1)
    C_j = a0 + Cr*exp(pi*1i*j/N);
    f_cj = f_eval(theta, mu, J_2, alpha, c, C_j, semi, r1, r2, delta_t);
    sum1 = sum1 + exp(2*pi*1i*j/N)*f_cj;
    sum2 = sum2 + exp(pi*1i*j/N)*f_cj;
end

f_plus = f_eval(theta, mu, J_2, alpha, c, a0 + Cr, semi, r1, r2, delta_t);
f_minus = f_eval(theta, mu, J_2, alpha, c, a0 - Cr, semi, r1, r2, delta_t);

numerator = f_plus + f_minus + 2*sum1;
denominator = f_plus - f_minus + 2*sum2;

a = a0 + Cr*real(numerator)/real(denominator);

%% Calculate Velocity
% Find the Velocity Magnitutde (Vis Viva)
v1_mag = sqrt(mu*(2/norm(r1) - 1/a));
v2_mag = sqrt(mu*(2/norm(r2) - 1/a));
% Find the orbital plane and time 1 and 2
[psi, ~] = get_angles(semi,c,a);
e = get_eccentricity(r1, r2, theta, psi, a);
[i, Omega_1, Omega_2] = newton_angles_2(r1, r2, a, e, J_2, mu, alpha, delta_t);


%% Results

% Absolute Error in velocity magnitude
v1_error = abs(v1_mag_real - v1_mag);
v2_error = abs(v2_mag_real - v2_mag);

% Absolute Error in RAAN
Omega1_error = abs(Omega_1 - Omega_1_v);
Omega2_error = abs(Omega_2 - Omega_2_v);

%% Functions
function [psi, cphi] = get_angles(semi,c,a)

    alpha = 2*asin(sqrt(semi/(2*a)));
    beta = 2*asin(sqrt((semi-c)/(2*a)));
    psi = (alpha - beta)/2;

    cphi = cos((alpha + beta)/2);

end

function e = get_eccentricity(r1_vec, r2_vec, theta, psi, a)
    r1 = norm(r1_vec);
    r2 = norm(r2_vec);
    
    % Calculate Orbit Parameter
    p = norm(r1)*norm(r2)*sin(theta/2)^2/(a*sin(psi)^2);

    % Claculate eccentricity
    e = sqrt(1 - p/a);
end

function [i, Omega_1, Omega_2] = newton_angles_2(r1_vec, r2_vec, a, e, J2, mu, alpha, delta_t)
    Beta = 3/2*((J2*sqrt(mu)*alpha^2)/(a^(7/2)*(1 - e^2)^2));
    
    % F = [r_1.h_hat, r2.h_hat, Omega_1 + delta_Omega - Omega_2]
    % X = [i; Omega_1, Omega_2]
    F = @(x) [sin(x(1))*sin(x(2))*r1_vec(1) - sin(x(1))*cos(x(2))*r1_vec(2) + cos(x(1))*r1_vec(3);
              sin(x(1))*sin(x(3))*r2_vec(1) - sin(x(1))*cos(x(3))*r2_vec(2) + cos(x(1))*r2_vec(3);
              x(3) + Beta*cos(x(1))*delta_t - x(2)];
    x = fsolve(F, [.5, 0, 0]);

    i = x(1);
    Omega_1 = x(2);
    Omega_2 = x(3);
end

function val = f_eval(theta, mu, J2, alpha, c, a, semi, r1, r2, delta_t)
    % find the angles
    [psi, cphi] = get_angles(semi, c, a);
    e = get_eccentricity(r1, r2, theta, psi, a);
    [i,~,~] = newton_angles_2(r1, r2, a, e, J2, mu, alpha, delta_t);
    val = (2*(psi - sin(psi)*cphi) + (J2*((alpha/(2*a))^2)*(3*sin(i)^2 - 2))/(((1 - e^2)^3))*(4*(e^2 + 2)*psi - 16*sin(psi)*cphi + 8*sin(2*psi)*cphi^2 - 4*e^2*sin(2*psi)) - sqrt(mu/(a^3))*delta_t)^-1;
end

