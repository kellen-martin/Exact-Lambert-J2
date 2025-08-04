function [psi, cphi] = lagrange_angles(semi,c,a, theta, ap_cross)
    alpha = 2*asin(sqrt(semi/(2*a)));
    beta = 2*asin(sqrt((semi-c)/(2*a)));
    if theta > pi
        beta = -beta;
    end
    if ap_cross
        alpha = 2*pi - alpha;
    end

    psi = (alpha - beta)/2;
    cphi = cos((alpha + beta)/2);

end