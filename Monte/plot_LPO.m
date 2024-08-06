clear all; close all; clc;

% plot_LPO.m
% Benjamin Hanson, 2024

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
x(:,1:2) = x(:,1:2)*const.LU;
x(:,3:4) = x(:,3:4)*const.LU/const.TU;

%% Plotting Nominal Trajectories
figure(1); 
plot(x(:,1),x(:,2), 'r-','LineWidth',1,'DisplayName','Nominal');
drawnow;
figure(2); 
plot(x(:,3),x(:,4), 'r-','LineWidth', 1, 'DisplayName','Nominal');
drawnow;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
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