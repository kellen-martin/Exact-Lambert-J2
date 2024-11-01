%% This Scipt is for testing the Generalized Labert's Problem 
clc;
clear
%% Given
% Position vectors
r1 = [8000; 0; 0];    % [km] 
r2 = [-7000; 8000; 500];    % [km]

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

%% Test Feval
test = f_eval(theta, mu, J_2, alpha, c, a0, semi, r1, r2, delta_t);
%% Trapaziod Rule

N = 6; % [number of points]
sum1 = 0;
sum2 = 0;
for j=1:(N-1)
    C_j = a0 + Cr*exp(pi*1i*j/N);
    sum1 = sum1 + exp(2*pi*1i*j/N)*f_eval(theta, mu, J_2, alpha, c, C_j, semi, r1, r2, delta_t);
    sum2 = sum2 + exp(pi*1i*j/N)*f_eval(theta, mu, J_2, alpha, c, C_j, semi, r1, r2, delta_t);
end

f_plus = f_eval(theta, mu, J_2, alpha, c, a0 + Cr, semi, r1, r2, delta_t);
f_minus = f_eval(theta, mu, J_2, alpha, c, a0 - Cr, semi, r1, r2, delta_t);

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

function [i, Omega_1, Omega_2] = newton_angles(r1_vec, r2_vec, a, e, J2, mu, alpha, delta_t)
    % Iterative Method for the inclination and RAAN
    % Initial Guess
    h_hat = cross(r1_vec,r2_vec)/norm(cross(r1_vec,r2_vec));
    i = asin(h_hat(3));
    f = @(Omega) sin(i)*cos(Omega)*r1_vec(1) + sin(i)*sin(Omega)*r1_vec(2) + cos(i)*r1_vec(3);
    % Omega_1 = fsolve(f,0);
    Omega_1 = pi;

    % Initialize Newton Method
    tol = 1E-10;
    n = 0;
    N_max = 100;
    F = [10; 10];

    % Newton Method
    while(norm(F)>tol && n<N_max)
        Beta = 3/2*((J2*sqrt(mu)*alpha^2)/(a^(7/2)*(1 - e^2)^2));
        delta_Omega = -Beta*cos(i)*delta_t;
        Omega_2 = Omega_1 + delta_Omega;
        ddeltaOmega_di = Beta*sin(i)*delta_t;

        F = [sin(i)*cos(Omega_1)*r1_vec(1) + sin(i)*sin(Omega_1)*r1_vec(2) + cos(i)*r1_vec(3);
            sin(i)*cos(Omega_2)*r2_vec(1) + sin(i)*sin(Omega_2)*r2_vec(2) + cos(i)*r2_vec(3)];

        DF = [cos(i)*sin(Omega_1)*r1_vec(1) + cos(i)*sin(Omega_1)*r1_vec(2) - sin(i)*r1_vec(3), -sin(i)*sin(Omega_1)*r1_vec(1) + sin(i)*cos(Omega_1)*r1_vec(2);
            (cos(i)*cos(Omega_2) - sin(i)*sin(Omega_2)*ddeltaOmega_di)*r2_vec(1) + (cos(i)*sin(Omega_2) + sin(i)*cos(Omega_2)*ddeltaOmega_di)*r2_vec(2) - sin(i)*r2_vec(3),-sin(i)*sin(Omega_2)*r2_vec(1) + sin(i)*cos(Omega_2)*r2_vec(2)];
        
        delta_X = DF\F;
        i = i + delta_X(1);
        Omega_1 = Omega_1 + delta_X(2); 
        n = n + 1;
    end
    if(n == N_max)
        fprintf('Did not Converge')
    end
end

function val = f_eval(theta, mu, J2, alpha, c, a, semi, r1, r2, delta_t)
    % find the angles
    [psi, cphi] = get_angles(semi, c, a);
    e = get_eccentricity(r1, r2, theta, psi, a);
    [i,~,~] = newton_angles(r1, r2, a, e, J2, mu, alpha, delta_t);
    val = sqrt(mu/(a^3))/(2*(psi - sin(psi)*cphi) + (J2*((alpha/(2*a))^2)*(3*sin(i)^2 - 2))/(((1 - e^2)^3))*(4*(e^2 + 2)*psi - 16*sin(psi)*cphi + 4*sin(2*psi)*(cphi^2 - e^2/2)));
end

