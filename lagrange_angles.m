function [psi, cphi] = lagrange_angles(semi,c,a)

    alpha = 2*asin(sqrt(semi/(2*a)));
    beta = 2*asin(sqrt((semi-c)/(2*a)));
    psi = (alpha - beta)/2;

    cphi = cos((alpha + beta)/2);

end