function v = ref_orbit(r, mu)
    rp_min = 6400;
    rp = 0;
    ra_max = 11000;
    ra = 1e9;
    n = 0;
    N_max = 10;
    while (rp<rp_min || ra>ra_max) && n < N_max
    v_dir = rand(3,1);
    v_hat = v_dir/norm(v_dir);

    a_max = 11000;

    e_max = .99;
    e = e_max*rand;

    a_min = 2*norm(r);
    a = a_min + (a_max - a_min)*rand;

    v_mag = sqrt(mu*(2/norm(r) - 1/a));
    v = v_mag*v_hat;
    h = cross(r,v);
    if h(3) < 0
        v = -v;
    end
    [a, e, ~, ~, ~, ~, ~] = get_oe(r, v, mu);
    rp = a*(1 - e);
    ra = a*(1-e);
    n = n+1;
    end

end