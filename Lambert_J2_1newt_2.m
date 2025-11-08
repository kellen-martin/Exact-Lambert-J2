%% J2 Perturbed Lambert Solver
function [a, v1, v2] = Lambert_J2_1newt_2(r1, r2, delta_t, mu, J_2, alpha, N)

%% Compute Geometric Quantities
% angle between position vectors
theta = acos(dot(r1/norm(r1), r2/norm(r2)));
check = cross(r1, r2);

if(check(3)<0)
    theta = 2*pi-theta;

end
% cord length
c = norm(r1 - r2);

% semi-perimeter
semi = .5*(norm(r1) + norm(r2) + c);

tminE = MinEtof(semi, c, theta);
ap_cross = delta_t>tminE;

%% Compute Contour
a_min_1 = semi/2;
e_max = .5;
r_min = min(norm(r1), norm(r2));
r_a = -r_min*(e_max + 1)/(e_max - 1);
a_max = (r_a + r_min)/2;

r_max = max(norm(r1), norm(r2));
r_p = -r_max*(e_max - 1)/(e_max + 1);
a_min_2 = (r_max + r_p)/2;
a_min = max(a_min_1, a_min_2);
Cr = (a_max - a_min)/2;
a0 = a_min + Cr;



n = 0;

%% Compute Inclination
% Find the Inclination for a = a0
for k=1:2
[psi, ~] = lagrange_angles(semi, c, a0, theta, ap_cross);
e = get_eccentricity(r1, r2, theta, psi, a0);
[i, ~, ~] = newton_angles(r1, r2, a0,e,J_2,mu,alpha,delta_t);

%% Composite Trapezoid
sum1 = 0;
sum2 = 0;
for j=1:(N-1)
    C_j = a0 + Cr*exp(pi*1i*j/N);
    f_cj = f_eval_1newt(theta, mu, J_2, alpha, c, C_j, semi, r1, r2, i, delta_t, ap_cross);
    sum1 = sum1 + exp(2*pi*1i*j/N)*f_cj;
    sum2 = sum2 + exp(pi*1i*j/N)*f_cj;
end

f_plus = f_eval_1newt(theta, mu, J_2, alpha, c, a0 + Cr, semi, r1, r2, i, delta_t, ap_cross);
f_minus = f_eval_1newt(theta, mu, J_2, alpha, c, a0 - Cr, semi, r1, r2, i,  delta_t, ap_cross);

numerator = f_plus + f_minus + 2*sum1;
denominator = f_plus - f_minus + 2*sum2;

a = a0 + Cr*real(numerator)/real(denominator);

a0 = a;
delta = 30;
a_min = max(a0 - delta, a_min);
a_max = a0 + delta;
Cr = (a_max - a_min)/2;
a0 = a_min + Cr;
end
%% Calculate Velocity
% Find the Velocity Magnitutde (Vis Viva)
v1_mag = sqrt(mu*(2/norm(r1) - 1/a));
v2_mag = sqrt(mu*(2/norm(r2) - 1/a));

% Find the orbital plane and time 1 and 2
[psi, ~] = lagrange_angles_v(semi,c,a, theta, ap_cross);
e = get_eccentricity(r1, r2, theta, psi, a);
p = get_parameter(r1, r2, theta, psi, a);
[i, Omega_1, Omega_2] = newton_angles(r1, r2, a, e, J_2, mu, alpha, delta_t);

[v1, v2] = velocity_solve(r1, r2, v1_mag, v2_mag, theta, Omega_1, Omega_2, i, p, e, mu);
end