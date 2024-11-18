%% J2 Perturbed Lambert Solver
function [a, v1, v2] = Lamabert_J2(r1, r2, delta_t, mu, J_2, alpha, N)

%% Compute Geometric Quantities
% angle between position vectors
theta = acos(dot(r1/norm(r1), r2/norm(r2)));

% cord length
c = norm(r1 - r2);

% semi-perimeter
semi = .5*(norm(r1) + norm(r2) + c);

%% Compute Contour
a_min = c/2;
e_max = .8;
r_p = min(norm(r1), norm(r2));
r_a = -r_p*(e_max + 1)/(e_max - 1);
a_max = (r_a + r_p)/2;
Cr = (a_max - a_min)/2;
a0 = a_min + Cr;

%% Composite Trapezoid
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
v1 = sqrt(mu*(2/norm(r1) - 1/a));
v2 = sqrt(mu*(2/norm(r2) - 1/a));

% % Find the orbital plane and time 1 and 2
% [psi, ~] = get_angles(semi,c,a);
% e = get_eccentricity(r1, r2, theta, psi, a);
% [i, Omega_1, Omega_2] = newton_angles_2(r1, r2, a, e, J_2, mu, alpha, delta_t);
end