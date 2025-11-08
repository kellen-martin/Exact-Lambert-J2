%% Test of Classical Orbit Element Integration using ode45
clc
clear
close all

mu = 3.986*10^5;    % [km^3/s^2]
J_2 = 1.0826E-3;
alpha = 6378;          % [km]
N = 16;

[r1_lc, v1_lc, ~, ~, delta_t, ~] = lambert_conditions(mu);


[a, e, p, i, Omega, omega, f] = get_oe(r1_lc, v1_lc, mu);
oes_0 = [a, e, i, omega, Omega, f];

delta_t = round(delta_t);
% delta_t = 24*60*60;
[oes, t1] = oeIntegrateJ2(oes_0, delta_t);


for i = 1:length(oes)
    [r(i, :), ~] = get_rv(oes(i,:), mu);
end
figure
plot3(r(:,1), r(:,2), r(:,3))
hold on
plot3(r(end,1), r(end,2), r(end,3), 'go')

[r1, v1] = get_rv(oes_0, mu);
[r2, v2] = get_rv(oes(end, :), mu);

[r2_L, v2_L, a_out] = lopezPropegate(r1, v1, delta_t);
[r2_V, ~] = pkepler(r1, v1, delta_t, 0, 0);

diff_L = norm(r2 - r2_L);
diff_V = norm(r2 - r2_V);

[a_s, v1_s, v2_s] = Lambert_J2_1newt_2(r1_lc, r2, delta_t, mu, J_2, alpha, N);

%% Check with integ
[a, e, p, inc, Omega, omega, f] = get_oe(r1, v1_s, mu);
oes_0s = [a, e, inc, omega, Omega, f];

[oes_s, ts] = oeIntegrateJ2(oes_0s, delta_t);

[r2_i_check, ~] = get_rv(oes_s(end, :), mu);
err = norm(r2_i_check - r2);

[r2_Ls, ~, ~] = lopezPropegate(r1, v1_s, delta_t);
err_L = norm(r2_Ls - r2);

%%
for i = 1:length(oes)
    [r(i, :), ~] = get_rv(oes(i,:), mu);
end
for i = 1:length(oes_s)
    [rs(i, :), ~] = get_rv(oes_s(i,:), mu);
end
figure
plot3(r(:,1), r(:,2), r(:,3))
hold on
plot3(r(end,1), r(end,2), r(end,3), 'go')
plot3(rs(:,1), rs(:, 2), rs(:, 3))
plot3(rs(end,1), rs(end, 2), rs(end, 3), 'r*')

c = norm(r1 - r2);
semi = .5*(norm(r1) + norm(r2) + c);
a_min=semi/2;
theta = acos(dot(r1/norm(r1), r2/norm(r2)));
check = cross(r1, r2);

if(check(3)<0)
    theta = 2*pi-theta;
end

tminE = MinEtof(semi, c, theta);
check = delta_t>tminE;


val = f_eval_1newt(theta, mu, J_2, alpha, c, a_s, semi, r1, r2, inc, delta_t, check)^-1;

%%
time = '+ ' + string(delta_t) + ' sec';

% Create an instance of STK
% Get reference to running STK instance
uiApplication = actxGetRunningServer('STK12.Application');

% Get our IAgStkObjectRoot interface
root = uiApplication.Personality2;
scen = root.CurrentScenario;

scen.SetTimePeriod('10 Jun 2025 00:00:00.000', time);

% Display Sats

% Shooting Method
J2_prop = 'ePropagatorHPOP';
HPOP = 'ePropagatorHPOP';

prop = J2_prop;

% Perturbed Solution Sat
    colr = 16711680;     % Blue
    sat = scen.Children.Item( 'Solution' );  % Get current satellite
    sat.SetPropagatorType( prop );  % already set in scenario
    sat.Propagator.Step = 10.0;
    sat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', r1(1), r1(2), r1(3), v1(1), v1(2), v1(3) );
    sat.Propagator.Propagate;
    
    % display satellite 
    graphics = sat.Graphics;
    graphics.SetAttributesType('eAttributesBasic');
    attributes = graphics.Attributes;
    attributes.Inherit = false;
    attributes.Line.Width = 3;
    %attributes.Line.Style = 'eLongDash';
    attributes.Color = colr;
    
    % Display one orbit pass and no ground track on 3D
    passdata3D = sat.VO.Pass.TrackData.PassData;
    groundTrack3D = passdata3D.GroundTrack;
    groundTrack3D.SetLeadDataType('eDataNone');
    groundTrack3D.SetTrailSameAsLead;
    orbit3D = passdata3D.Orbit;
    orbit3D.SetLeadDataType('eDataOnePass');
    orbit3D.SetTrailSameAsLead;

% Error
satPosDP = sat.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(scen.StartTime,scen.StopTime,60);
satx = cell2mat(satPosDP.DataSets.GetDataSetByName('x').GetValues);
saty = cell2mat(satPosDP.DataSets.GetDataSetByName('y').GetValues);
satz = cell2mat(satPosDP.DataSets.GetDataSetByName('z').GetValues);
r2_check = [satx(end); saty(end); satz(end)];

error_bad = r2 - r2_check;
abs_error_test = norm(error_bad);