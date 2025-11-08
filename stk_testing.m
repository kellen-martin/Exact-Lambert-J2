%% Find the error
clc;
clear

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

N = 256;

%% Set Initial Conditions
trials = 100;

r1_save = zeros(3,trials);
r2_save = zeros(3,trials);
delta_t_save = zeros(1,trials);

h = waitbar(0, 'Running');
for i=1:trials
[r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);
r1_save(:,i) = r1;
r2_save(:,i) = r2;
delta_t_save(:,i) = delta_t;

time = '+ ' + string(delta_t) + ' sec';
%% Lambert Solvers
    [~, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);
    % tol = 1e-9;
    % [v_1v,v_2v, nv] = UV_Lambert_Bisect(r1, r2, mu, delta_t, tol, false);


%% Create an instance of STK
% Get reference to running STK instance
uiApplication = actxGetRunningServer('STK12.Application');

% Get our IAgStkObjectRoot interface
root = uiApplication.Personality2;
scen = root.CurrentScenario;

scen.SetTimePeriod('10 Jun 2025 00:00:00.000', time);

%% Display Sats

% Shooting Method
J2_prop = 'ePropagatorJ2Perturbation';
HPOP = 'ePropagatorHPOP';

prop = J2_prop;

% Perturbed Solution Sat
    colr = 16711680;     % Blue
    sat = scen.Children.Item( 'Solution' );  % Get current satellite
    sat.SetPropagatorType( prop );  % already set in scenario
    sat.Propagator.Step = 10.0;
    sat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', r1(1), r1(2), r1(3), v1_L(1), v1_L(2), v1_L(3) );
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

% Unpert, solution sat

    % colr = 16777215;     %  white
    % unsat = scen.Children.Item( 'un' );  % Get current satellite
    % unsat.SetPropagatorType( prop );  % already set in scenario
    % unsat.Propagator.Step = 10.0;
    % unsat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', r1(1), r1(2), r1(3), v_1v(1), v_1v(2), v_1v(3) );
    % unsat.Propagator.Propagate;
    % 
    % % display satellite 
    % graphics = unsat.Graphics;
    % graphics.SetAttributesType('eAttributesBasic');
    % attributes = graphics.Attributes;
    % attributes.Inherit = false;
    % attributes.Line.Width = 3;
    % %attributes.Line.Style = 'eLongDash';
    % attributes.Color = colr;
    % 
    % % Display one orbit pass and no ground track on 3D
    % passdata3D = unsat.VO.Pass.TrackData.PassData;
    % groundTrack3D = passdata3D.GroundTrack;
    % groundTrack3D.SetLeadDataType('eDataNone');
    % groundTrack3D.SetTrailSameAsLead;
    % orbit3D = passdata3D.Orbit;
    % orbit3D.SetLeadDataType('eDataOnePass');
    % orbit3D.SetTrailSameAsLead;



%% Error
satPosDP = sat.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(scen.StartTime,scen.StopTime,60);
satx = cell2mat(satPosDP.DataSets.GetDataSetByName('x').GetValues);
saty = cell2mat(satPosDP.DataSets.GetDataSetByName('y').GetValues);
satz = cell2mat(satPosDP.DataSets.GetDataSetByName('z').GetValues);
r2_check = [satx(end); saty(end); satz(end)];

error = r2 - r2_check;
abs_error(i) = norm(error);

% satPosDP = unsat.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(scen.StartTime,scen.StopTime,60);
% unsatx = cell2mat(satPosDP.DataSets.GetDataSetByName('x').GetValues);
% unsaty = cell2mat(satPosDP.DataSets.GetDataSetByName('y').GetValues);
% unsatz = cell2mat(satPosDP.DataSets.GetDataSetByName('z').GetValues);
% r2_checkun = [unsatx(end); unsaty(end); unsatz(end)];
% 
% unerror = r2 - r2_checkun;
% abs_errorun(i) = norm(unerror);

waitbar(i/trials, h)
end

%% Post Process
abs_err_mean = mean(abs_error);
% abs_err_meanun = mean(abs_errorun);
figure
histogram(abs_error)

[abs_err_rm, rm_idx] = rmoutliers(abs_error,"grubbs");
figure
histogram(abs_err_rm);


c = vecnorm(r1_save - r2_save, 2);
semi = .5*(vecnorm(r1_save,2) + vecnorm(r2_save,2) + c);
a_min=semi/2;

for i=1:trials
    r1i = r1_save(:,i);
    r2i = r2_save(:,i);
    theta(i) = acos(dot(r1i/norm(r1i), r2i/norm(r2i)));
    check = cross(r1i, r2i);

    if(check(3)<0)
        theta(i) = 2*pi-theta(i);
    end
end

delta_t_rm = delta_t_save(~rm_idx);
theta_rm = theta(~rm_idx);
a_min_rm = a_min(~rm_idx);

figure
scatter(delta_t_rm, abs_err_rm, '*')
title('Time-of-Flight vs Error')
xlabel('Time-of-Flight [sec]')
ylabel('Absolute Error [km]')

figure
scatter(theta_rm, abs_err_rm, '*')
title('Time-of-Flight vs Error')
xlabel('Time-of-Flight [sec]')
ylabel('Absolute Error [km]')

%% get a's
for i=1:trials
    r1i = r1_save(:,i);
    r2i = r2_save(:,i);
    delta_ti = delta_t_save(i);

[a(i), ~, ~] = Lambert_J2_1newt_2(r1i, r2i, delta_ti, mu, J_2, alpha, N);
end

%% Plots
a_min_close = (a - a_min)./a;

figure
scatter(a_min_close, abs_error, '*')
title('Error near a_{min}')
xlabel('Releative Difference between Solution and a_{min}')
ylabel('Absolute Error [km]')

idx = find(a_min_close > .02);
err_rm = abs_error(idx);
err_rm_mean = mean(err_rm);
tof_rm = delta_t_save(idx);
theta_rm = theta(idx);

figure
histogram(err_rm, 100, 'Normalization','probability')
hold on
xline(err_rm_mean)
hold off
title('Absolute Error Distribution')
xlabel('Absolute Error')

figure
scatter(tof_rm, err_rm, '*')
title('Time-of-Flight vs Error')
xlabel('Time-of-Flight [sec]')
ylabel('Absolute Error [km]')

figure
scatter(theta_rm, err_rm, '*')
title('Transfer Angle vs Error')
xlabel('Transfer Angle [rad]')
ylabel('Absolute Error [km]')

figure
scatter(a_min_close(idx), err_rm, '*')
title('Error near a_{min}')
xlabel('Releative Difference between Solution and a_{min}')
ylabel('Absolute Error [km]')

%% Remove bad thetas
idx1 = 2.8 < theta_rm;
idx2 = theta_rm < 3.2;
idx3 = theta_rm>5.74;

idxt = ~((idx1 & idx2) | idx3);

err_rmt = err_rm(idxt);
err_rm_meant = mean(err_rmt);
tof_rmt = tof_rm(idxt);
theta_rmt = theta_rm(idxt);

figure
histogram(err_rmt, 100, 'Normalization','probability')
hold on
xline(err_rm_meant)
hold off
title('Absolute Error Distribution')
xlabel('Absolute Error')

figure
scatter(tof_rmt, err_rmt, '*')
title('Time-of-Flight vs Error')
xlabel('Time-of-Flight [sec]')
ylabel('Absolute Error [km]')

figure
scatter(theta_rmt, err_rmt, '*')
title('Transfer Angle vs Error')
xlabel('Transfer Angle [rad]')
ylabel('Absolute Error [km]')


% %% More Things
% [bad, bi] = max(abs_error);
% 
% r1_bad = r1_save(:,bi);
% r2_bad = r2_save(:,bi);
% delta_t_bad = delta_t_save(bi);
% 
% time = '+ ' + string(delta_t) + ' sec';
% % Lambert Solvers
%     [~, v1_L, v2_L] = Lamabert_J2_1newt(r1_bad, r2_bad, delta_t_bad, mu, J_2, alpha, N);
%     % tol = 1e-9;
%     % [v_1v,v_2v, nv] = UV_Lambert_Bisect(r1, r2, mu, delta_t, tol, false);
% 
% 
% % Create an instance of STK
% % Get reference to running STK instance
% uiApplication = actxGetRunningServer('STK12.Application');
% 
% % Get our IAgStkObjectRoot interface
% root = uiApplication.Personality2;
% scen = root.CurrentScenario;
% 
% scen.SetTimePeriod('10 Jun 2025 00:00:00.000', time);
% 
% % Display Sats
% 
% % Shooting Method
% J2_prop = 'ePropagatorJ2Perturbation';
% HPOP = 'ePropagatorHPOP';
% 
% prop = J2_prop;
% 
% % Perturbed Solution Sat
%     colr = 16711680;     % Blue
%     sat = scen.Children.Item( 'Solution' );  % Get current satellite
%     sat.SetPropagatorType( prop );  % already set in scenario
%     sat.Propagator.Step = 10.0;
%     sat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', r1(1), r1(2), r1(3), v1_L(1), v1_L(2), v1_L(3) );
%     sat.Propagator.Propagate;
% 
%     % display satellite 
%     graphics = sat.Graphics;
%     graphics.SetAttributesType('eAttributesBasic');
%     attributes = graphics.Attributes;
%     attributes.Inherit = false;
%     attributes.Line.Width = 3;
%     %attributes.Line.Style = 'eLongDash';
%     attributes.Color = colr;
% 
%     % Display one orbit pass and no ground track on 3D
%     passdata3D = sat.VO.Pass.TrackData.PassData;
%     groundTrack3D = passdata3D.GroundTrack;
%     groundTrack3D.SetLeadDataType('eDataNone');
%     groundTrack3D.SetTrailSameAsLead;
%     orbit3D = passdata3D.Orbit;
%     orbit3D.SetLeadDataType('eDataOnePass');
%     orbit3D.SetTrailSameAsLead;
% 
% % Unpert, solution sat
% 
%     % colr = 16777215;     %  white
%     % unsat = scen.Children.Item( 'un' );  % Get current satellite
%     % unsat.SetPropagatorType( prop );  % already set in scenario
%     % unsat.Propagator.Step = 10.0;
%     % unsat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', r1(1), r1(2), r1(3), v_1v(1), v_1v(2), v_1v(3) );
%     % unsat.Propagator.Propagate;
%     % 
%     % % display satellite 
%     % graphics = unsat.Graphics;
%     % graphics.SetAttributesType('eAttributesBasic');
%     % attributes = graphics.Attributes;
%     % attributes.Inherit = false;
%     % attributes.Line.Width = 3;
%     % %attributes.Line.Style = 'eLongDash';
%     % attributes.Color = colr;
%     % 
%     % % Display one orbit pass and no ground track on 3D
%     % passdata3D = unsat.VO.Pass.TrackData.PassData;
%     % groundTrack3D = passdata3D.GroundTrack;
%     % groundTrack3D.SetLeadDataType('eDataNone');
%     % groundTrack3D.SetTrailSameAsLead;
%     % orbit3D = passdata3D.Orbit;
%     % orbit3D.SetLeadDataType('eDataOnePass');
%     % orbit3D.SetTrailSameAsLead;
% 
% 
% 
% % Error
% satPosDP = sat.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(scen.StartTime,scen.StopTime,60);
% satx = cell2mat(satPosDP.DataSets.GetDataSetByName('x').GetValues);
% saty = cell2mat(satPosDP.DataSets.GetDataSetByName('y').GetValues);
% satz = cell2mat(satPosDP.DataSets.GetDataSetByName('z').GetValues);
% r2_check = [satx(end); saty(end); satz(end)];
% 
% error_bad = r2 - r2_check;
% abs_error_bad = norm(error_bad);
% 
% % satPosDP = unsat.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(scen.StartTime,scen.StopTime,60);
% % unsatx = cell2mat(satPosDP.DataSets.GetDataSetByName('x').GetValues);
% % unsaty = cell2mat(satPosDP.DataSets.GetDataSetByName('y').GetValues);
% % unsatz = cell2mat(satPosDP.DataSets.GetDataSetByName('z').GetValues);
% % r2_checkun = [unsatx(end); unsaty(end); unsatz(end)];
% % 
% % unerror = r2 - r2_checkun;
% % abs_errorun(i) = norm(unerror);
% 
