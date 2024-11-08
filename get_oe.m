function [a, e, p, i, Omega, omega, f] = get_oe(r, rdot, mu)
% ECI unit vectors
k_hat = [0 0 1];
i_hat = [1 0 0];
j_hat = [0 1 0];

% Angular momentum
h = cross(r, rdot);
h_hat = h/norm(h);

% Eccentricity
e_vec = (1/mu)*cross(rdot, h) - r/norm(r);
e_hat = e_vec/norm(e_vec);
e_hat_perp = cross(h_hat, e_hat);
e = norm(e_vec);

% Orbit parameter
p = norm(h)^2/mu;

% Semi-Major axis
a = p/(1 - e^2);

% Inclination
i = acos(dot(h_hat, k_hat));

% Node Vector
n_hat = cross(k_hat, h_hat)/norm(cross(k_hat, h_hat));
n_hat_perp = cross(h_hat, n_hat);

% Longitude of Ascending node
Omega = atan2(dot(j_hat, n_hat),dot(i_hat, n_hat));

% Arguement of periapsis
omega = atan2(dot(e_vec, n_hat_perp), dot(e_vec, n_hat));

% True Anomaly
f = atan2(dot(r, e_hat_perp), dot(r, e_hat));


end 