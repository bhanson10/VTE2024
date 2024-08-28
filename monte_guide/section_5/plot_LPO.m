clear all; close all; clc;

% plot_LPO.m
% Benjamin Hanson, 2024

%% Initial Conditions
const.mu = 2.528017528540000E-5; const.LU = 668519; const.TU = 48562;  
const.U = [const.LU, const.LU, const.LU/const.TU, const.LU/const.TU];
const.T = 2.6513344042156235E+0; const.r =  1560.8*500;
rv.start = [1.017714765; -1.069793E-20; -1.197784E-13; 1.187104E-2];

%% Analytical
tspan = [0, const.T];
x0 = rv.start;
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
[t, x] = ode87(@(t, x) PCR3BP(t, x, const), tspan, x0, options);

ax(:,1:4) = x(:,1:4).*const.U;
t = t.*const.TU;

%% Initializing Figures
lbl.XString = '$x$ (km)'; lbl.YString = '$y$ (km)';
initialize_figures('n', 1, 'spacing', {[50 100 700 700]},'lbl', lbl, 'axs', {'equal'}); 
set(gca, 'FontName', 'Times', 'FontSize', 14);

% Loading Europa
const.mu = 2.52801752854E-5; const.LU = 668519; const.TU = 48562; const.cV = [const.LU, const.LU, const.LU/const.TU, const.LU/const.TU]; 
europa_2d = imread('./Europa_2D.png');
europa_2d = imresize(europa_2d, 1.5);
[height, width, ~] = size(europa_2d);
image([(1-const.mu)*const.LU-width/2, (1-const.mu)*const.LU+width/2], [-height/2, height/2], europa_2d);

lbl.XString = '$v_x$ (km/s)';
lbl.YString = '$v_y$ (km/s)';
initialize_figures('n', 2, 'spacing', {[800 100 700 700]}, 'lbl', lbl, 'axs', {'equal'}); 
set(gca, 'FontName', 'Times', 'FontSize', 14);

%% Truth
[x, n] = parse_monte_txt("jupiter_europa_lpo.txt");
mx(:,1:4) = x(:,1:4).*const.U;

%% Plotting Analytical Trajectories
figure(1); 
plot(ax(:,1),ax(:,2), 'b:','LineWidth',2,'DisplayName','Analytical');
drawnow;
figure(2); 
plot(ax(:,3),ax(:,4), 'b:','LineWidth', 2, 'DisplayName','Analytical');
drawnow;

%% Plotting Monte Trajectories
figure(1);
plot(mx(:,1),mx(:,2),'r--','LineWidth',2,'DisplayName','Monte');
legend()
drawnow;
figure(2);
plot(mx(:,3),mx(:,4),'r--','LineWidth', 2, 'DisplayName','Monte');
legend()
drawnow;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = PCR3BP(t, x, const)
    x1 = [x(3); x(4); 2*x(4)+x(1)-(const.mu*(x(1)-1+const.mu)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*(x(1)+const.mu)/(((x(1)+const.mu)^2+x(2)^2)^(1.5))); -2*x(3)+x(2)-(const.mu*x(2)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*x(2)/(((x(1)+const.mu)^2+x(2)^2)^(1.5)))];
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, n] = parse_monte_txt(filename)
    fileID = fopen(filename, 'r');
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        x(count, :) = [str2double(line{1}); str2double(line{2}); str2double(line{4}); str2double(line{5})];
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(x); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%