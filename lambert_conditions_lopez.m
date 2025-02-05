function [r1, v1, r2, v2, delta_t] = lambert_conditions_lopez()
%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

% Time step
t_step = 1;        % [s]

%% Random Initial Contitions
[oes, delta_t] = random_orbit(mu);
a = oes(1);
e = oes(2);
i = oes(3);
omega = oes(4);
Omega = oes(5);

% Compute Mean anomaly from true anomaly
f = oes(6);
M = atan2(-sqrt(1 - e^2)*sin(f), -e - cos(f)) + pi - e*(sqrt(1 - e^2)*sin(f))/(1 + e*cos(f));
data_oes = [a e i omega Omega M];

%% Write datos.problema
folder_path = 'C:\Users\kmartin6\Desktop\ppkbj21-Kellen\ppkbj21-Kellen';
return_path = 'C:\Users\kmartin6\Desktop\Research\Exact-Lambert-J2';
cd(folder_path);

fileID_data = fopen('datos.problema', 'w');

fprintf(fileID_data, 'O\n\n');
fprintf(fileID_data, '%.16g ', data_oes);    % Orbital Elements [a e i omega Omega M]
fprintf(fileID_data, '\n\n');
fprintf(fileID_data, '0.0\n\n');             % Initial Time
fprintf(fileID_data, '%.16g', delta_t);      % Final Time
fprintf(fileID_data, '\n\n');
fprintf(fileID_data, '%.16g', t_step);       % Time step
fprintf(fileID_data, '\n\n');
fprintf(fileID_data, '%.16g', mu);           % Gravitational Parameter
fprintf(fileID_data, '\n\n');
fprintf(fileID_data, '%.16g', alpha);        % Earth Mean Radius
fprintf(fileID_data, '\n\n');
fprintf(fileID_data, '0.0\n\n'); 
fprintf(fileID_data, '%.16g', J_2);          % Earth J_2 constant

%% Compile and Execute ppkb
files = dir(fullfile(folder_path,'*.c'));
files_list = strjoin(fullfile(folder_path, {files.name}, ''));
compile_cmd = ['gcc ' files_list ' -O3 -std=c89 -lm -o ppkb'];

system(compile_cmd);
system('ppkb.exe');

%% Read first and last entry of solcart.out
fileID_oe_out = fopen('solorb.out','r');
fileID_cart_out = fopen('solcart.out','r');

oe_out = fscanf(fileID_oe_out, '%f %f %f %f %f %f %f \n', [7 Inf]);
cart_out = fscanf(fileID_cart_out, '%f %f %f %f %f %f %f \n', [7 Inf]);

r1 = cart_out(2:4,1);
r2 = cart_out(2:4,end);
v1 = cart_out(5:7,1);
v2 = cart_out(5:7,end);


% return to original folder
cd(return_path)
end