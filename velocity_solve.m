function [v1, v2] = velocity_solve(r1, r2, v1_mag, v2_mag, Omega_1, Omega_2, inc, p, mu)
% Find specific angular momentum and 
h_mag = sqrt(p*mu);
h1_hat = [sin(inc)*sin(Omega_1); -sin(inc)*cos(Omega_1); cos(inc)];
h2_hat = [sin(inc)*sin(Omega_2); -sin(inc)*cos(Omega_2); cos(inc)];

% Find flight path angle
gama_1 = asin(h_mag/(norm(r1)*v1_mag));
gama_2 = asin(h_mag/(norm(r2)*v2_mag));

% Find basis for the orbital planes
r1_hat = r1/norm(r1);
r2_hat = r2/norm(r2);

theta1_hat = cross(h1_hat, r1_hat);
theta2_hat = cross(h2_hat, r2_hat);

% Find the velocity vector
v1 = v1_mag*(cos(gama_1)*r1_hat + sin(gama_1)*theta1_hat);
v2 = v2_mag*(cos(gama_2)*r2_hat + sin(gama_2)*theta2_hat);
end