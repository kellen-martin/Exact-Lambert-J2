function [r, v, oes, energies, hs] = integrate_2bp_j2(r0, v0, mu, time)
State0 = [r0; v0];

% Givens 
R_0 = 6400; % [km]
J_2 = 1E-3; 

% ODE Call 
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
tspan = [0 time];
[t,State] = ode45(@(t, State) ode_fun_J2(State, mu), tspan, State0, options);

% Decompose State Vector
r = State(:, 1:3);
v = State(:, 4:6);

% Compute orbital elements, energy, and angular momentum
oes = zeros(length(v),7);
energies = zeros(1,length(v));
hs = zeros(1,length(v));
for j = 1:length(v)
[a, e, p, i, Omega, omega, f] = get_oe(r(j,:), v(j,:), mu);
oes(j,:) = [a, e, p, i, Omega, omega, f];
energies(j) = .5*norm(v(j,:)).^2 - mu/norm(r(j,:)) - mu/(2*norm(r(j,:))^3)*R_0.^2*J_2*(1 - 3*r(j,3).^2/norm(r(j,:))^2);
hs(j) = dot(cross(r(j,:), v(j,:)),[0; 0; 1]);
end
end