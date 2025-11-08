%% Test of Solving
clc
clear
close all

mu = 3.986*10^5;    % [km^3/s^2]
J_2 = 1.0826E-3;
alpha = 6378;          % [km]
Ns = [8, 16, 32, 64, 128, 256];
% Ns = [32];

trials = 1000;

h = waitbar(0, 'Running');
for n = 1:length(Ns)
    N = Ns(n);
    for k = 1:trials
        isBad = true;
        while isBad
        [r1_lc, v1_lc, ~, ~, delta_t, ~] = lambert_conditions(mu);


        [a, e, p, i, Omega, omega, f] = get_oe(r1_lc, v1_lc, mu);
        oes_0 = [a, e, i, omega, Omega, f];

        delta_t = round(delta_t);
        oes = oeIntegrateJ2(oes_0, delta_t);


        [r1, v1] = get_rv(oes_0, mu);
        [r2, v2] = get_rv(oes(end, :), mu);


        [a_s, v1_s, v2_s] = Lambert_J2_1newt_2(r1_lc, r2, delta_t, mu, J_2, alpha, N);

        %% Check with integ
        c = norm(r1 - r2);
        semi = .5*(norm(r1) + norm(r2) + c);
        a_min=semi/2;
        theta(k) = acos(dot(r1/norm(r1), r2/norm(r2)));
        check = cross(r1, r2);

        if(check(3)<0)
            theta(k) = 2*pi-theta(k);
        end

        a_min_close(k) = (a_s - a_min)/a_s;
        if a_min_close(k) < .02 || a_min>a_s || (theta(k)>2.9 && theta(k)<3.2) || theta(k) > 5.8
        else
        isBad = false;
        [a, e, p, inc, Omega, omega, f] = get_oe(r1, v1_s, mu);
        % oes_0s = [a, e, inc, omega, Omega, f];
        % 
        % oes_s = oeIntegrateJ2(oes_0s, delta_t);
        % 
        % [r2_i_check, ~] = get_rv(oes_s(end, :), mu);
        % err = norm(r2_i_check - r2);
        end
        end
    tminE = MinEtof(semi, c, theta(k));
    check = delta_t>tminE;


    val_1(k, n) = abs(f_eval_1newt(theta(k), mu, J_2, alpha, c, a_s, semi, r1, r2, inc, delta_t, check)^-1);
    end
    waitbar((k + trials*(n-1))/(trials*length(Ns)), h)

    a_mean(n) = mean(val_1(:, n));
end
close(h)

%% Process Data
bin_edges = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
figure
histogram(val_1(:,1), bin_edges)
hold on
set(gca, 'XScale', 'log')
for i = 2:length(Ns)
    histogram(val_1(:, i), bin_edges)
end
hold off 
val_means = mean(val_1, 1);

figure
plot(Ns, val_means, 'k*-')
xlabel('N (Quaderature Points)')
ylabel('Absolute Error')
title('Equation Solution Convergence')

figure
scatter(a_min_close,val_1)

figure
scatter(theta, val_1)