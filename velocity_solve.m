function [v1, v2] = velocity_solve(r1, r2, v1_mag, v2_mag, Omega_1, Omega_2, inc, theta)
h_mag = v1_mag*norm(r1)*sin(theta);
h_1 = h_mag*[sin(inc)*sin(Omega_1), sin(inc)*cos(Omega_1),cos(inc)];
h_2 = h_mag*[sin(inc)*sin(Omega_2), sin(inc)*cos(Omega_2),cos(inc)];

% We want to solve the system r x v = h
F1 = @(v1) [h_1(1) - (r1(2)*v1_mag*v1(3) - r1(3)*v1_mag*v1(2));
            h_1(2) + (r1(1)*v1_mag*v1(3) - r1(3)*v1_mag*v1(1));
            h_1(3) - (r1(1)*v1_mag*v1(2) - r1(2)*v1_mag*v1(1))];

F2 = @(v2) [h_2(1) - (r2(2)*v2_mag*v2(3) - r2(3)*v2_mag*v2(2));
            h_2(2) + (r2(1)*v2_mag*v2(3) - r2(3)*v2_mag*v2(1));
            h_2(3) - (r2(1)*v2_mag*v2(2) - r2(2)*v2_mag*v2(1))];

v1 = v1_mag*fsolve(F1, [0, 0, 0]);
v2 = v2_mag*fsolve(F2, [0, 0, 0]);
end