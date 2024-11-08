%% This Scipt is for testing the Generalized Labert's Problem 
clc;
clear
%% Given

% time of flight
delta_t = 5000;        % [s]   

% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

% Position vectors
r1 = [7000; 0; 0];    % [km] 

%% Test Case
% Initial Velocity
v1 = [0; 7.5; 5];

% Propegation 
%[r2,v2] = pkepler(r1, v1, delta_t, ndot, nddot);
[r_prop, v_prop] = integrate_2bp_j2(r1, v1, mu, delta_t);

% Final State
r2 = r_prop(end,:);
v2 = v_prop(end,:);

% Calculate Velocity Magnitudes
v1_mag_real = norm(v1);
v2_mag_real = norm(v2);

% Inital and Final RAAN and inclination
[a_prop, e_prop, ~, i_prop, Omega_1_prop, ~, ~] = get_oe(r_prop(1,:), v1, mu);
[~, ~, ~, ~, Omega_2_prop, ~, ~] = get_oe(r_prop(end,:), v2, mu);

% Plot Test Case
figure
plot3(r_prop(:,1), r_prop(:,2), r_prop(:,3))
hold on
axis equal
[x,y,z] = sphere;
x = alpha*x;
y = alpha*y;
z = alpha*z;
surf(x,y,z)
hold off

%% Test thing
[i, Omega_1, Omega_2] = newton_angles_2(r1, r2, a_prop, e_prop, J_2, mu, alpha, delta_t)
%% Compute other quantites
% angle between r1 & r2
theta = dot(r1, r2)/(norm(r1)*norm(r2));

% cord length
c = norm(r1 - r2);

% semi perimeter
semi = .5*(norm(r1) + norm(r2) + c);

%% Compute contour
a_min = 1000;
e_max = .9;
r_p = min(norm(r1), norm(r2));
r_a = -r_p*(e_max + 1)/(e_max - 1);
a_max = (r_a + r_p)/2;
Cr = (a_max - a_min)/2;
a0 = a_min + Cr;


%% Trapaziod Rule

N = 100; % [number of points]
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

a = a0 + Cr*real((f_plus + f_minus + 2*sum1)/(f_plus + f_minus + 2*sum2));

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
Omega1_error = abs(Omega_1 - Omega_1_prop);
Omega2_error = abs(Omega_2 - Omega_2_prop);
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
function [i, Omega_1, Omega_2] = newton_angles_2(r1_vec, r2_vec, a, e, J2, mu, alpha, delta_t)
    Beta = 3/2*((J2*sqrt(mu)*alpha^2)/(a^(7/2)*(1 - e^2)^2));
    
    % F = [r_1.h_hat, r2.h_hat, Omega_1 + delta_Omega - Omega_2]
    % X = [i; Omega_1, Omega_2]
    F = @(x) [sin(x(1))*cos(x(2))*r1_vec(1) + sin(x(1))*sin(x(2))*r1_vec(2) + cos(x(1))*r1_vec(3);
              sin(x(1))*cos(x(3))*r2_vec(1) + sin(x(1))*sin(x(3))*r2_vec(2) + cos(x(1))*r2_vec(3);
              x(2) - Beta*cos(x(1))*delta_t - x(3)];
    x = fsolve(F, [.5, 0, 0]);

    i = x(1);
    Omega_1 = x(2);
    Omega_2 = x(3);
end

function [i, Omega_1, Omega_2] = newton_angles(r1_vec, r2_vec, a, e, J2, mu, alpha, delta_t)
    % Iterative Method for the inclination and RAAN
    % Initial Guess

    i = asin(r1_vec(3)/norm(r1_vec));
    f = @(Omega) sin(i)*cos(Omega)*r1_vec(1) + sin(i)*sin(Omega)*r1_vec(2) + cos(i)*r1_vec(3);
    Omega_1 = fsolve(f,0);


    % Initialize Newton Method
    tol = 1E-10;
    n = 0;
    N_max = 1000;
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
        fprintf('Did not Converge \n')
    end
end

function val = f_eval(theta, mu, J2, alpha, c, a, semi, r1, r2, delta_t)
    % find the angles
    [psi, cphi] = get_angles(semi, c, a);
    e = get_eccentricity(r1, r2, theta, psi, a);
    [i,~,~] = newton_angles_2(r1, r2, a, e, J2, mu, alpha, delta_t);
    val = (2*(psi - sin(psi)*cphi) + (J2*((alpha/(2*a))^2)*(3*sin(i)^2 - 2))/(((1 - e^2)^3))*(4*(e^2 + 2)*psi - 16*sin(psi)*cphi + 4*sin(2*psi)*(cphi^2 - e^2/2)) - sqrt(mu/(a^3))*delta_t);
end

