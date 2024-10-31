%% This Scipt is for testing the Generalized Labert's Problem 
clc;
clear
%% Given
% Position vectors
r1 = [8000, 0, 0];    % [km] 
r2 = [-7000, 8000, 500];    % [km]

% time of flight
delta_t = 900000;        % [s]   

% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

%% Compute other quantites
% angle between r1 & r2
theta = dot(r1, r2)/(norm(r1)*norm(r2));

% direction of the angular momentum vecotr
h_hat = cross(r1,r2)/norm(cross(r1,r2));

% s = sin^2(i)
s = h_hat(3)^2;

% cord length
c = norm(r1 - r2);

% semi perimeter
semi = .5*(norm(r1) + norm(r2) + c);

%% Compute contour
a_min = alpha;
e_max = .9;
r_p = min(norm(r1), norm(r2));
r_a = -r_p*(e_max + 1)/(e_max - 1);
a_max = (r_a + r_p)/2;
Cr = (a_max - a_min)/2;
a0 = a_min + Cr;

%% Trapaziod Rule

N = 6; % [number of points]
sum1 = 0;
sum2 = 0;
for j=1:(N-1)
    C_j = a0 + Cr*exp(pi*1i*j/N);
    sum1 = sum1 + exp(2*pi*1i*j/N)*f_eval(theta, mu, J_2, alpha, s, c, C_j, semi, r1, r2);
    sum2 = sum2 + exp(pi*1i*j/N)*f_eval(theta, mu, J_2, alpha, s, c, C_j, semi, r1, r2);
end

f_plus = f_eval(theta, mu, J_2, alpha, s, c, a0 + Cr, semi, r1, r2);
f_minus = f_eval(theta, mu, J_2, alpha, s, c, a0 - Cr, semi, r1, r2);

a = a0 + Cr*real(f_plus + f_minus + 2*sum1)/real(f_plus + f_minus + 2*sum2);

%% Calculate Velocity

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

    e = sqrt(1 - (r1*r2*sin(theta/2)^2)/(a^2*sin(psi)^2));
end

function [i, Omega] = newton_angles(r1_vec, r2_vec, a, e, J2, mu, delta_t)
    % Iterative Method for the inclination and RAAN

    i = 0;
    Omega = 0;
end

function val = f_eval(theta, mu, J2, alpha, s, c, a, semi, r1, r2)
    % find the angles
    [psi, cphi] = get_angles(semi, c, a);
    e = get_eccentricity(r1, r2, theta, psi, a);
    val = sqrt(mu/(a^3))/(2*(psi - sin(psi)*cphi) + (J2*((alpha/(2*a))^2)*(3*s - 2))/(((1 - e^2)^3))*(4*(e^2 + 2)*psi - 16*sin(psi)*cphi + 4*sin(2*psi)*(cphi^2 - e^2/2)));
end

