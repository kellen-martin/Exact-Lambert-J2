%% Trajectory optimization with J2
% clc;
% clear
% close all

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

tol = 1e-6;
max_iter = 1e6;
trials = 5e5;
diff = zeros([trials, 1]);
dt_out = zeros([trials, 1]);
r1_out = zeros([3, trials]);
v1_out = zeros([3, trials]);
r2_out = zeros([3, trials]);
v2_out = zeros([3, trials]);
rng('shuffle')
h = waitbar(0, 'Running');
%% Generate Random Orbits
for i = 1:trials
    psi = -1;
    while psi<0
        [r1, v1_sol, r2, v2_sol, dt1, oes] = lambert_conditions(mu);
        v1 = ref_orbit(r1, mu);
        v2 = ref_orbit(r2, mu);
        [v1p, v2m, psi] = UV_Lambert2(r1, r2, mu, dt1, tol, max_iter, false);
    end
    delta_v1 = v1p - v1;
    delta_v2 = v2 - v2m;
    delta_v = norm(delta_v1) + norm(delta_v2);

    N = 256;
    [a, v1p_j, v2m_j] = Lamabert_J2_1newt(r1, r2, dt1, mu, J_2, alpha, N);
    delta_v1_j = v1p_j - v1;
    delta_v2_j = v2 - v2m_j;
    delta_v_j = norm(delta_v1_j) + norm(delta_v2_j);

    if isreal(v1p_j)
    [rc,vc] = pkepler(r1, v1p_j, dt1, 0, 0);
    check = norm(r2 - rc)/norm(r2);
    else
        check = 1e9;
    end

    if check > .03
        %fprintf("bad solution")
        diff(i) = nan;
    else
        diff(i) = delta_v_j - delta_v;
    end
    dt_out(i) = dt1;
    r1_out(:,i) = r1;
    v1_out(:,i) = v1;
    r2_out(:,i) = r2;
    v2_out(:,i) = v2;
    waitbar(i/trials, h)
end
close(h)
[largest_diff, index] = min(diff);
r1 = r1_out(:, index);
v1 = v1_out(:, index);
r2 = r2_out(:, index);
v2 = v2_out(:, index);
dt = dt_out(index);
[v1p, v2m, psi] = UV_Lambert2(r1, r2, mu, dt, tol, max_iter, false);
[a, v1p_j, v2m_j] = Lamabert_J2_1newt(r1, r2, dt, mu, J_2, alpha, N);

ref1_a = get_oe(r1, v1, mu);
t1 = 2*pi*sqrt(ref1_a^3/mu);
t1_vec = linspace(0, t1, 500);
ref2_a = get_oe(r2, v2, mu);
t2 = 2*pi*sqrt(ref2_a^3/mu);
t2_vec = linspace(0, t2, 500);
t_transfer_vec = linspace(0, dt, 500);

num_improved = length(find(diff<-1));

for i = 1:length(t_transfer_vec)
    [ref1(:,i), ~] = pkepler(r1, v1, t1_vec(i), 0, 0);
    [ref2(:,i), ~] = pkepler(r2, v2, t2_vec(i), 0, 0);
    [rn(:,i), ~] = pkepler(r1, v1p, t_transfer_vec(i), 0, 0);
    [rj(:,i), ~] = pkepler(r1, v1p_j, t_transfer_vec(i), 0, 0);
end

% Make Earth
r = alpha;
[x,y,z] = sphere(50);
x = x*r;
y = y*r;
z = z*r;

% figure
plot3(rn(1,:), rn(2,:), rn(3,:))
hold on
plot3(rj(1,:), rj(2,:), rj(3,:))
plot3(ref1(1,:), ref1(2,:), ref1(3,:))
plot3(ref2(1,:), ref2(2,:), ref2(3,:))
plot3(r1(1), r1(2), r1(3), 'b*')
plot3(r2(1), r2(2), r2(3), 'r*')
surf(x,y,z, 'FaceColor',"#77AC30", 'EdgeColor','none', 'FaceAlpha',.5)
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
legend('Classical Solution', 'Perturbed Solution', 'Departure Orbit', 'Arrival Orbit')

%% Post Processing
mu_d = nanmean(diff);
hist_data = diff(~isnan(diff));
sigma_d = std(hist_data);

g1 = find(hist_data<=-1);
g5 = find(hist_data<=-5);
g10 = find(hist_data<=-10);
g15 = find(hist_data<=-15);
g20 = find(hist_data<=-20);

numg1 = length(g1);
numg5 = length(g5);
numg10 = length(g10);
numg15 = length(g15);
numg20 = length(g20);

%% Trajectory plots
g1l5 = find(diff >= 15 & diff <= 20);
g5l10 = find(diff <= -5 & diff > -10);
g10l15 = find(diff <= -10 & diff > -15);
g15l20 = find(diff <= -15 & diff > -20);

index_1 = g1l5(1);
diff_1 = diff(index_1);

index_5 = g5l10(1);
diff_5 = diff(index_5);

index_10 = g10l15(1);
diff_10 = diff(index_10);

index_15 = g15l20(1);
diff_15 = diff(index_15);

%% Plot at proper index
% input proper index
r1 = r1_out(:, index_1);
v1 = v1_out(:, index_1);
r2 = r2_out(:, index_1);
v2 = v2_out(:, index_1);
dt = dt_out(index_1);
theta = acos(dot(r1,r2)/(norm(r1*norm(r2))))
[v1p, v2m, psi] = UV_Lambert2(r1, r2, mu, dt, tol, max_iter, false);
[a, v1p_j, v2m_j] = Lamabert_J2_1newt(r1, r2, dt, mu, J_2, alpha, N);

ref1_a = get_oe(r1, v1, mu);
t1 = 2*pi*sqrt(ref1_a^3/mu);
t1_vec = linspace(0, t1, 500);
ref2_a = get_oe(r2, v2, mu);
t2 = 2*pi*sqrt(ref2_a^3/mu);
t2_vec = linspace(0, t2, 500);
t_transfer_vec = linspace(0, dt, 500);

num_improved = length(find(diff<-1));

for i = 1:length(t_transfer_vec)
    [ref1(:,i), ~] = pkepler(r1, v1, t1_vec(i), 0, 0);
    [ref2(:,i), ~] = pkepler(r2, v2, t2_vec(i), 0, 0);
    [rn(:,i), ~] = pkepler(r1, v1p, t_transfer_vec(i), 0, 0);
    [rj(:,i), ~] = pkepler(r1, v1p_j, t_transfer_vec(i), 0, 0);
end

% Make Earth
r = alpha;
[x,y,z] = sphere(50);
x = x*r;
y = y*r;
z = z*r;

figure
plot3(rn(1,:), rn(2,:), rn(3,:))
hold on
plot3(rj(1,:), rj(2,:), rj(3,:))
plot3(ref1(1,:), ref1(2,:), ref1(3,:))
plot3(ref2(1,:), ref2(2,:), ref2(3,:))
plot3(r1(1), r1(2), r1(3), 'b*')
plot3(r2(1), r2(2), r2(3), 'r*')
surf(x,y,z, 'FaceColor',"#77AC30", 'EdgeColor','none', 'FaceAlpha',.5)
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
legend('Classical Solution', 'Perturbed Solution', 'Departure Orbit', 'Arrival Orbit')