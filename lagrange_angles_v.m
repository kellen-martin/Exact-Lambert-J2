function [psi, cphi] = lagrange_angles_v(semi,c,a, theta)
% Lagrange angles for calculating the velocity
% remove imaginary part of alpha if close to a_min
    alpha = 2*asin(sqrt(semi/(2*a)));
    alpha = real(alpha);
    beta = 2*asin(sqrt((semi-c)/(2*a)));

    psi = (alpha - beta)/2;

    cphi = cos((alpha + beta)/2);

end