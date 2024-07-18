clear all; close all; clc;

% plot_gbees.m
% Benjamin Hanson, 2024

%% Define colors
load("colors.mat"); 

%% Initial Condition
traj.sys = '/PCR3BP'; name = 'LPO_1';
filename = append('./ICs', traj.sys,'/', name,'.csv');
po = readmatrix(filename); const.mu = po(8);
traj.rv0 = [po(2); po(3); po(5); po(6)]; 
const.LU = 668519; const.TU = 48562; const.mu = po(12); const.cV = [const.LU, const.LU, const.LU/const.TU, const.LU/const.TU]; 
const.T = po(9); const.r = po(10)*500;
colors = [0 0 0.1724; 1 0 0; 0 0 1; 0 1 0]; 

%% Initializing Figures
lbl.XString = '$x$ (km)'; lbl.YString = '$y$ (km)';
initialize_figures('n', 1, 'spacing', {[50 100 700 700]},'lbl', lbl, 'axs', {'equal'}); 
set(gca, 'FontName', 'Times', 'FontSize', 14);

% Loading Europa
europa_2d = imread('./Europa_2D.png');
europa_2d = imresize(europa_2d, 1.5);
[height, width, ~] = size(europa_2d);
image([(1-const.mu)*const.LU-width/2, (1-const.mu)*const.LU+width/2], [-height/2, height/2], europa_2d);

lbl.XString = '$v_x$ (km/s)';
lbl.YString = '$v_y$ (km/s)';
initialize_figures('n', 2, 'spacing', {[800 100 700 700]}, 'lbl', lbl, 'axs', {'equal'}); 
set(gca, 'FontName', 'Times', 'FontSize', 14);

%% Truth
num_dist = 16; 
measure_t = const.T/num_dist; 
tspan = [0:measure_t:const.T]; 
t = 0; x0 = traj.rv0; x = x0'; t_idx = 1;
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
count = 1; 
for i=tspan(2:end)
    [ti, xi] = ode87(@(t, x) f(t, x, const), [tspan(count),i], x0, options);
    x = [x; xi(2:end,:)]; t = [t; ti(2:end)]; t_idx = [t_idx; numel(t)];
    x0 = x(end,:)'; count = count + 1; 
end

x(:,1:4) = x(:,1:4).*const.cV; 
t = t.*const.TU; 

%% Plotting Nominal Trajectories
figure(1); 
plot(x((t <= 14*3600),1),x((t <= 14*3600),2), 'Color', tolred, 'LineStyle', '-','LineWidth',2,'DisplayName','Nominal');
plot(x(:,1),x(:,2), 'Color', tolred, 'LineStyle', '--','LineWidth',1,'DisplayName','Nominal');
drawnow;
figure(2); 
plot(x((t <= 14*3600),3),x((t <= 14*3600),4), 'Color', tolred, 'LineStyle', '-','LineWidth',2,'DisplayName','Nominal');
plot(x(:,3),x(:,4), 'Color', tolred, 'LineStyle', '--','LineWidth', 1, 'DisplayName','Nominal');
drawnow;

%{
%% Particle Filter
dir_path = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/PF/cmake-build-debug/Epochs/Europa/M0";
file_list = dir(fullfile(dir_path, '*.txt'));  % List only .txt files
num_files = numel(file_list);

count = 1; 
for i=[0, num_files-1]
    file_path = dir_path + "/pdf_" + num2str(i) + ".txt"; 

    [x_pf, P_pf, n_pf, t_pf(count)] = parse_nongaussian_txt(file_path);

    % Converting units
    x_pf(:,1:4) = x_pf(:,1:4).*const.cV; 
    t_pf(count) = t_pf(count)*const.TU;

    xest_pf{count} = zeros(size(x_pf(1,:)));
    for j=1:n_pf
        xest_pf{count} = xest_pf{count}+x_pf(j,:).*P_pf(j);
    end
    
    figure(1); 
    scatterHandle = scatter(x_pf(:,1), x_pf(:,2), 9, [0.8, 0.8 0.8], 'filled');
    scatterHandle.ZData = -1 * ones(size(x_pf(:,1))); 
    figure(2); 
    scatterHandle = scatter(x_pf(:,3), x_pf(:,4), 9, [0.8, 0.8 0.8], 'filled');
    scatterHandle.ZData = -1 * ones(size(x_pf(:,3))); 
    drawnow; 

    count = count + 1;
end
%}

%% GBEES
NM = 1; 
clear p; p.color = tolblue; p.name = "GBEES"; p.alpha = [0.2, 0.4, 0.6]; 
dir_path = "./gbees/v0/Data/PDFs/P";

count = 1;
for nm=0:NM-1

    sub_dir_path = dir_path + num2str(nm); 
    file_list = dir(fullfile(sub_dir_path, '*.txt'));  % List only .txt files
    num_files = numel(file_list);

    for i=[0, num_files-1]
        file_path = sub_dir_path + "/pdf_" + num2str(i) + ".txt"; 
    
        [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(file_path);
    
        % Converting units
        x_gbees(:,1:4) = x_gbees(:,1:4).*const.cV; 
        t_gbees(count) = t_gbees(count)*const.TU;
    
        xest_gbees{count} = zeros(size(x_gbees(1,:)));
        for j=1:n_gbees
            xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
        end

        figure(1); 
        plot_nongaussian_surface(x_gbees(:,1:2),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
        figure(2); 
        plot_nongaussian_surface(x_gbees(:,3:4),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
        drawnow; 

        count = count + 1;
    end
end

%{
%% Comparison Parameters
% Position/Velocity Error
n = size(xest_pf); n = n(2);
pe = NaN(n,3); ve = NaN(n,3);
for i=1:n
    xi_pf = xest_pf{i}; 
    xi_ukf = xest_ukf{i};
    xi_enkf = xest_enkf{i};
    xi_gbees = xest_gbees{i};

    pe(i,1) = norm((x(t_idx(i),1:2)'-xi_ukf(1:2)));
    ve(i,1) = norm((x(t_idx(i),3:4)'-xi_ukf(3:4)));
    pe(i,2) = norm((x(t_idx(i),1:2)'-xi_enkf(1:2)));
    ve(i,2) = norm((x(t_idx(i),3:4)'-xi_enkf(3:4)));
    pe(i,3) = norm((x(t_idx(i),1:2)-xi_gbees(1:2)));
    ve(i,3) = norm((x(t_idx(i),3:4)-xi_gbees(3:4)));
end

% Jaccard Coefficient
n = size(xest_pf); n = n(2); 
po_ukf = NaN(n,3); vo_ukf = NaN(n,3); 
po_enkf = NaN(n,3); vo_enkf = NaN(n,3);
po_gbees = NaN(n,3); vo_gbees = NaN(n,3);
for i=1:n
    shpi_pf = shps_pf{i}; 
    shpi_ukf = shps_ukf{i}; 
    shpi_enkf = shps_enkf{i}; 
    shpi_gbees = shps_gbees{i}; 

    meani_ukf = xest_ukf{i}; covi_ukf = Sest_ukf{i};
    meani_enkf = xest_enkf{i}; covi_enkf = Sest_enkf{i};
    
    pi_pf = shpi_pf{1}; vi_pf = shpi_pf{2}; 
    pi_ukf = shpi_ukf{1}; vi_ukf = shpi_ukf{2}; 
    pi_enkf = shpi_enkf{1}; vi_enkf = shpi_enkf{2}; 
    pi_gbees = shpi_gbees{1}; vi_gbees = shpi_gbees{2}; 

    for j=1:numel(pi_ukf)
        p_ukf_shp.sd = sd(j); 

        p_ukf_shp.Points = pi_ukf{j}'; 
        p_ukf_shp.mean = meani_ukf(1:2); 
        p_ukf_shp.cov = covi_ukf(1:2,1:2); 
        po_ukf(i,j) = jaccard_coef(pi_pf{j}, p_ukf_shp);

        v_ukf_shp.sd = sd(j); 

        v_ukf_shp.Points = vi_ukf{j}'; 
        v_ukf_shp.mean = meani_ukf(3:4); 
        v_ukf_shp.cov = covi_ukf(3:4,3:4); 
        vo_ukf(i,j) = jaccard_coef(vi_pf{j}, v_ukf_shp);
    end

    for j=1:numel(pi_enkf)
        p_enkf_shp.sd = sd(j); 

        p_enkf_shp.Points = pi_enkf{j}'; 
        p_enkf_shp.mean = meani_enkf(1:2); 
        p_enkf_shp.cov = covi_enkf(1:2,1:2); 
        po_enkf(i,j) = jaccard_coef(pi_pf{j}, p_enkf_shp);

        v_enkf_shp.sd = sd(j); 

        v_enkf_shp.Points = vi_enkf{j}'; 
        v_enkf_shp.mean = meani_enkf(3:4); 
        v_enkf_shp.cov = covi_enkf(3:4,3:4); 
        vo_enkf(i,j) = jaccard_coef(vi_pf{j}, v_enkf_shp);
    end
    
    for j=1:numel(pi_gbees)
        po_gbees(i,j) = jaccard_coef(pi_pf{j}, pi_gbees{4-j});
        vo_gbees(i,j) = jaccard_coef(vi_pf{j}, vi_gbees{4-j});
    end
end

% Timing
file_path = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/PF/cmake-build-debug/Epochs/Europa/timing_1.txt";
[pt_pf, st_pf] = parse_timing_txt(file_path);
st_pf = st_pf.*const.TU; 

file_path = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/UKF/cmake-build-debug/Epochs/Europa/timing.txt";
[pt_ukf, st_ukf] = parse_timing_txt(file_path);
st_ukf = st_ukf.*const.TU; 

file_path = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/EnKF/cmake-build-debug/Epochs/Europa/timing.txt";
[pt_enkf, st_enkf] = parse_timing_txt(file_path);
st_enkf = st_enkf.*const.TU; 

file_path = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/GBEES/cmake-build-debug/Epochs/Europa/timing.txt";
[pt_gbees, st_gbees] = parse_timing_txt(file_path);
st_gbees = st_gbees.*const.TU;

pt_ukf = pt_ukf./pt_pf; pt_ukf(1) = 1;  
pt_enkf = pt_enkf./pt_pf; pt_enkf(1) = 1; 
pt_gbees = pt_gbees./pt_pf; pt_gbees(1) = 1; 

fnew = figure; hold on; fnew.Position = [50 150 1400 800]; 
tiledlayout(3, 1, "TileSpacing",'compact');

ax1 = nexttile; hold on;legend('Location','northwest','Interpreter','latex');
plot(t_pf,pe(:,1),'-o', 'Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF')
plot(t_pf,pe(:,2),'-o', 'Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF')
plot(t_pf,pe(:,3),'-o', 'Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','GBEES')

ax2 = nexttile; hold on;legend('Location','southwest','Interpreter','latex');
plot(t_pf,po_ukf(:,1),'-o','Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF - 1$\sigma$')
plot(t_pf,po_ukf(:,2),'-square','Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF - 2$\sigma$')
plot(t_pf,po_ukf(:,3),'-diamond','Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF - 3$\sigma$')
plot(t_pf,po_enkf(:,1),'-o','Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF - 1$\sigma$')
plot(t_pf,po_enkf(:,2),'-square','Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF - 2$\sigma$')
plot(t_pf,po_enkf(:,3),'-diamond','Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF - 3$\sigma$')
plot(t_pf,po_gbees(:,1),'-o','Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','GBEES - 1$\sigma$')
plot(t_pf,po_gbees(:,2),'-square','Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','GBEES - 2$\sigma$')
plot(t_pf,po_gbees(:,3),'-diamond','Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','GBEES - 3$\sigma$')

ax3 = nexttile; hold on;legend('Location','northwest','Interpreter','latex');
plot(st_ukf,pt_ukf,'-o', 'Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF')
plot(st_enkf,pt_enkf,'-o', 'Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF')
plot(st_gbees,pt_gbees,'-o', 'Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','GBEES')

tickvals = [0:0.5:round(t_pf(end)/3600)];
ticklabels = {}; count = 1;  for i=tickvals, ticklabels{count} = num2str(i); count = count + 1; end

set(ax1 , ...
     'FontName' , 'Times' , ...
     'FontSize' , 10 , ... 
     'XLim', [t_pf(1),1.001*t_pf(end)], ...
     'TickDir' , 'out' , ...
     'TickLength' , [0.0125 0.0125] , ...
     'XTick', [0:0.5*3600:round(t_pf(end)/3600)*3600], ...
     'XTickLabel', ticklabels, ...
     'XGrid' , 'on' , 'YGrid' , 'on' )  
ylabel( ax1 , "Position Error (km)", ...
                'Interpreter', 'LaTex' , 'Rotation', 90 ); 
ax1.YLabel.HorizontalAlignment = 'center';          
set(ax2 , ...
     'FontName' , 'Times' , ...
     'FontSize' , 10 , ... 
     'XLim', [t_pf(1),1.001*t_pf(end)], ...
     'TickDir' , 'out' , ...
     'TickLength' , [0.0125 0.0125] , ...
     'XTick', [0:0.5*3600:round(t_pf(end)/3600)*3600], ...
     'XTickLabel', ticklabels, ...
     'XGrid' , 'on' , 'YGrid' , 'on' ) 
ylabel(ax2, "$J_{\mathbf{r}}$", ...
                'Interpreter', 'LaTex' , 'Rotation', 90 ); 
ax2.YLabel.HorizontalAlignment = 'center';    
set(ax3 , ...
     'FontName' , 'Times' , ...
     'FontSize' , 10 , ... 
     'XLim', [t_pf(1),1.001*t_pf(end)], ...
     'YLim', [0,inf], ...
     'TickDir' , 'out' , ...
     'TickLength' , [0.0125 0.0125] , ...
     'XTick', [0:0.5*3600:round(t_pf(end)/3600)*3600], ...
     'XTickLabel', ticklabels, ...
     'XGrid' , 'on' , 'YGrid' , 'on' ) 
xlabel(ax3, "Time Past Epoch (hours)", ...
                'Interpreter', 'LaTex' , 'Rotation', 0 ); 
ylabel(ax3, "Normalized computation time", ...
                'Interpreter', 'LaTex' , 'Rotation', 90 ); 
ax3.YLabel.HorizontalAlignment = 'center'; 
%}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(t, x, const)
    x1 = [x(3); x(4); 2*x(4)+x(1)-(const.mu*(x(1)-1+const.mu)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*(x(1)+const.mu)/(((x(1)+const.mu)^2+x(2)^2)^(1.5))); -2*x(3)+x(2)-(const.mu*x(2)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*x(2)/(((x(1)+const.mu)^2+x(2)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        x(count, :) = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5})];
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, S, t] = parse_gaussian_txt(filename)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID)); 
    
    fgetl(fileID); % skip blank line

    line = split(fgetl(fileID)); % Read a line as a string
    x = [str2double(line{1});str2double(line{2});str2double(line{3});str2double(line{4})];
    
    fgetl(fileID); % skip blank line

    for i=1:4
        line = split(fgetl(fileID)); % Read a line as a string
        S(i,:) = [str2double(line{1});str2double(line{2});str2double(line{3});str2double(line{4})];
    end

    % Close the file
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pt, st] = parse_timing_txt(file_path)
    fileID = fopen(file_path, 'r');
    pt = []; st = [];
    
    while ~feof(fileID)
        line = split(fgetl(fileID)); 
        pt(end+1) = str2double(line{3});
        st(end+1) = str2double(line{7});
    end
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%