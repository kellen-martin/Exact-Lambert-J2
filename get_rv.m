function [r, rdot] = get_rv(oe, mu)
% Define orbit elements
a = oe(1);
e = oe(2);
i = oe(3);
omega = oe(4);
Omega = oe(5);
f = oe(6);
p = a*(1 - e^2);

% Perifocal Frame 
rmag_p = p/(1+e*cosd(f));
r_p = [rmag_p*cosd(f), rmag_p*sind(f), 0]';
rdot_p = sqrt(mu/p)*[-sind(f), e+cosd(f), 0]';

% Rotate to ECI 
C_np = angle2dcm(Omega, i, omega, 'ZXZ'); % DCM from N to P 
C_pn = transpose(C_np);                % DCM from P to N
r = mtimes(C_pn, r_p);
rdot = mtimes(C_pn, rdot_p);
end