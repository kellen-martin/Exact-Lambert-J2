function tof_minE = MinEtof(s, c, theta)
    mu = 3.986*10^5;
    alpha = pi;
    a = s/2;

    dm = 1;
    if theta > pi
        dm = -1;
    end

    beta =dm*2*asin(sqrt((s - c)/s));

    tof_minE = sqrt(a^3/mu)*(alpha - dm*(beta - sin(beta)));
end