clear all; close all; clc;

% pcr3bp_filter_analysis.m
% Benjamin Hanson, 2024

%% Initial Condition
rng(2024);
traj.sys = 'PCR3BP'; name = '1';
filename = append('./ICs/', traj.sys,'/IC_', name,'.csv');
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
imageFile = '/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/Europa_Map.png';
europaImg = imread(imageFile);
[X, Y, Z] = sphere(100);  
X = X.*const.r + (1-const.mu)*const.LU; Y = Y.*const.r; Z = Z.*const.r;
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', europaImg, 'EdgeColor', 'none', 'HandleVisibility', 'off');
scatter(1.02046139*const.LU,0,75,'d','k','LineWidth', 2, 'HandleVisibility', 'off');
drawnow;

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

x(:,1:2) = x(:,1:2).*const.LU; 
x(:,3:4) = x(:,3:4).*(const.TU/const.LU); 
t = t.*const.TU; 

%% Plotting Nominal Trajectories
figure(1); 
plot(x((t <= 14*3600),1),x((t <= 14*3600),2),'r-','LineWidth',2,'DisplayName','Nominal');
plot(x(:,1),x(:,2),'r--','LineWidth',1,'DisplayName','Nominal');
%scatter3(x(1,1),x(2,1),x(3,1),50,'filled','MarkerFaceColor','r','HandleVisibility','off');
drawnow;
figure(2); 
plot(x((t <= 14*3600),3),x((t <= 14*3600),4),'r-','LineWidth',2,'DisplayName','Nominal');
plot(x(:,3),x(:,4),'r--','LineWidth', 1, 'DisplayName','Nominal');
%scatter3(x(4,1),x(5,1),x(6,1),50,'filled','MarkerFaceColor','r','HandleVisibility','off');
drawnow;

%% Particle Filter
p.type = 'sharp'; p.color = colors(1,:); p.plt = 0; p.name = "PF"; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/PF/cmake-build-debug/Epochs/Europa/M1";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; 
for i=[0:numFiles-1]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x_pf, P_pf, n_pf, t_pf(count)] = parse_nongaussian_txt(FILE_PATH);

    % Converting units
    x_pf(:,1:4) = x_pf(:,1:4).*const.cV; 
    t_pf(count) = t_pf(count)*const.TU;

    xest_pf{count} = zeros(size(x_pf(1,:)));
    for j=1:n_pf
        xest_pf{count} = xest_pf{count}+x_pf(j,:).*P_pf(j);
    end
    
    [surfs_p] = plot_nongaussian_surface(x_pf(:,1:2),P_pf,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
    [surfs_v] = plot_nongaussian_surface(x_pf(:,3:4),P_pf,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
    surfs_pf{count} = {surfs_p, surfs_v}; 
    
    %{
    if(mod(i,4)==0)
        p.plt = 1; 
        axs{axs_count} = plot_non_gaussian(x_pf,P_pf,t_pf(count),p,axs_count,[0.30,0.92,0.99925],0.3,[]);
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end
    %}

    count = count + 1;
end

%% UKF
clear p; p.color = colors(2,:); p.plt = 0; p.name = "UKF";
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/UKF/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; sd = [1,2,3];
for i=[0:numFiles-1]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt";

    [xesti_ukf, Sesti_ukf, ti_ukf] = parse_gaussian_txt(FILE_PATH);

    xesti_ukf(:) = xesti_ukf(:).*const.cV(:); 
    Sesti_ukf    = Sesti_ukf.*(const.cV'*const.cV); 
    ti_ukf       = ti_ukf*const.TU; 

    xest_ukf{count} = xesti_ukf; Sest_ukf{count} = Sesti_ukf; t_ukf(count) = ti_ukf;
    
    [surfs_p] = plot_gaussian_ellipsoid(xesti_ukf(1:2,:),Sesti_ukf(1:2,1:2),sd,p);
    [surfs_v] = plot_gaussian_ellipsoid(xesti_ukf(3:4,:),Sesti_ukf(3:4,3:4),sd,p);
    surfs_ukf{count} = {surfs_p, surfs_v}; 

    %{
    if(mod(i,4)==0)
        p.plt = 1; 
        plot_gaussian(xest_ukf{count},Sest_ukf{count},t_ukf(count),p,axs_count,sd,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end
    %}

    count = count + 1;
end

%% EnKF
clear p; p.color = colors(3,:); p.plt = 0; p.name = "EnKF";
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/EnKF/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; sd = [1,2,3];
for i=[0:numFiles-1]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt";

    [xesti_enkf, Sesti_enkf, ti_enkf] = parse_gaussian_txt(FILE_PATH);

    xesti_enkf(:) = xesti_enkf(:).*const.cV(:); 
    Sesti_enkf    = Sesti_enkf.*(const.cV'*const.cV); 
    ti_enkf       = ti_enkf*const.TU; 
    
    xest_enkf{count} = xesti_enkf; Sest_enkf{count} = Sesti_enkf; t_enkf(count) = ti_enkf;

    [surfs_p] = plot_gaussian_ellipsoid(xesti_enkf(1:2,:),Sesti_enkf(1:2,1:2),sd,p);
    [surfs_v] = plot_gaussian_ellipsoid(xesti_enkf(3:4,:),Sesti_enkf(3:4,3:4),sd,p);
    surfs_enkf{count} = {surfs_p, surfs_v}; 
    
    %{
    if(mod(i,4)==0)
        p.plt = 1; 
        plot_gaussian(xest_enkf{count},Sest_enkf{count},t_enkf(count),p,axs_count,sd,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end
    %}

    count = count + 1;
end

%% GBEES
clear p; p.color = colors(4,:); p.plt = 0; p.name = "GBEES"; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/GBEES/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1;
for i=[0:numFiles-1]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(FILE_PATH);

    % Converting units
    x_gbees(:,1:4) = x_gbees(:,1:4).*const.cV; 
    t_gbees(count) = t_gbees(count)*const.TU;

    xest_gbees{count} = zeros(size(x_gbees(1,:)));
    for j=1:n_gbees
        xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
    end
    
    [surfs_p] = plot_nongaussian_surface(x_gbees(:,1:2),P_pf,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
    [surfs_v] = plot_nongaussian_surface(x_gbees(:,3:4),P_pf,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
    % [~,shp_p] = plot_grid_surface(x_gbees(:,1:2),P_gbees,[0.011, 0.125, 0.58],0.3,p);
    % [~,shp_v] = plot_grid_surface(x_gbees(:,3:4),P_gbees,[0.011, 0.125, 0.58],0.3,p);
    surfs_gbees{count} = {surfs_p, surfs_v}; 
    
    %{
    if(mod(i,4)==0)
        p.plt = 1; 
        plot_gbees(x_gbees,P_gbees,t_gbees(count),p,axs_count,[0.011, 0.125, 0.58],0.3,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end
    %}

    count = count + 1;
end

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
FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/PF/cmake-build-debug/Epochs/Europa/timing_1.txt";
[pt_pf, st_pf] = parse_timing_txt(FILE_PATH);
st_pf = st_pf.*const.TU; 

FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/UKF/cmake-build-debug/Epochs/Europa/timing.txt";
[pt_ukf, st_ukf] = parse_timing_txt(FILE_PATH);
st_ukf = st_ukf.*const.TU; 

FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/EnKF/cmake-build-debug/Epochs/Europa/timing.txt";
[pt_enkf, st_enkf] = parse_timing_txt(FILE_PATH);
st_enkf = st_enkf.*const.TU; 

FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/GBEES/cmake-build-debug/Epochs/Europa/timing.txt";
[pt_gbees, st_gbees] = parse_timing_txt(FILE_PATH);
st_gbees = st_gbees.*const.TU;

pt_ukf = pt_ukf./pt_pf; pt_ukf(1) = 1;  
pt_enkf = pt_enkf./pt_pf; pt_enkf(1) = 1; 
pt_gbees = pt_gbees./pt_pf; pt_gbees(1) = 1; 
%}

fnew = figure; hold on; fnew.Position = [50 150 1400 800]; 
tiledlayout(3, 1, "TileSpacing",'compact');

ax1 = nexttile; hold on;legend('Location','northwest','Interpreter','latex');
plot(t_pf,pe(:,1),'-o', 'Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF')
plot(t_pf,pe(:,2),'-o', 'Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF')
plot(t_pf,pe(:,3),'-o', 'Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','GBEES')

%{
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
%}

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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(t, x, const)
    x1 = [x(3); x(4); 2*x(4)+x(1)-(const.mu*(x(1)-1+const.mu)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*(x(1)+const.mu)/(((x(1)+const.mu)^2+x(2)^2)^(1.5))); -2*x(3)+x(2)-(const.mu*x(2)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*x(2)/(((x(1)+const.mu)^2+x(2)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axs] = plot_gbees(x,P,t,p,count,isovalues,alpha,axs)
    if isempty(axs)
        p.display = 1; 

        fi = figure(2*count+1); hold on; fi.Position = [50 100 700 700]; legend; 
        set(gca, 'FontName' , 'Times','FontSize',12); 
        plot_grid_surface(x(:,1:2),P,isovalues,alpha,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        title(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize', 18, 'FontName', 'Times', 'Interpreter','latex');
        axs_p = gca;
    
        fj = figure(2*count+2); hold on; fj.Position = [800 100 700 700]; legend; 
        set(gca, 'FontName' , 'Times','FontSize',12); 
        plot_grid_surface(x(:,3:4),P,isovalues,alpha,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        title(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize', 18, 'FontName', 'Times', 'Interpreter','latex');
        axs_v = gca; 

        axs = {axs_p,axs_v};
    else
        axs_p = axs{1}; axs_v = axs{2}; 
        
        p.display = 1; 
        
        p.axh = axs_p; 
        plot_grid_surface(x(:,1:2),P,isovalues,alpha,p);
        
        p.axh = axs_v; 
        plot_grid_surface(x(:,3:4),P,isovalues,alpha,p);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axs] = plot_non_gaussian(x,P,t,p,count,isovalues,alpha,axs)
    if isempty(axs)
        p.display = 1; 

        fi = figure(2*count+1); hold on; fi.Position = [50 100 700 700]; legend; 
        set(gca, 'FontName' , 'Times','FontSize',12); 
        plot_nongaussian_surface(x(:,1:2),P,isovalues,alpha,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        title(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize', 18, 'FontName', 'Times', 'Interpreter','latex');
        axs_p = gca;

        fj = figure(2*count+2); hold on; fj.Position = [800 100 700 700]; legend; 
        set(gca, 'FontName' , 'Times','FontSize',12); 
        plot_nongaussian_surface(x(:,3:4),P,isovalues,alpha,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        title(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize', 18, 'FontName', 'Times', 'Interpreter','latex');
        axs_v = gca; 

        axs = {axs_p,axs_v};
    else
        axs_p = axs{1}; axs_v = axs{2}; 
        
        p.display = 1; 
        
        p.axh = axs_p; 
        plot_nongaussian_surface(x(:,1:2),P,isovalues,alpha,p);
        
        p.axh = axs_v; 
        plot_nongaussian_surface(x(:,3:4),P,isovalues,alpha,p);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axs] = plot_gaussian(x,S,t,p,count,sd,axs)
    if isempty(axs)
        p.display=1; 

        fi = figure(2*count+1); hold on; fi.Position = [50 100 700 700]; legend; 
        set(gca, 'FontName' , 'Times','FontSize',12); 
        plot_gaussian_ellipsoid(x(1:2),S(1:2,1:2),sd,p);   
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        title(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize', 18, 'FontName', 'Times', 'Interpreter','latex');
        axs_p = gca;
    

        fj = figure(2*count+2); hold on; fj.Position = [800 100 700 700]; legend; 
        set(gca, 'FontName' , 'Times','FontSize',12); 
        plot_gaussian_ellipsoid(x(3:4),S(3:4,3:4),sd,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        title(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize', 18, 'FontName', 'Times', 'Interpreter','latex');
        axs_v = gca;

        axs = {axs_p,axs_v};
    else
        axs_p = axs{1}; axs_v = axs{2};  

        p.display = 1;

        p.axh = axs_p; 
        plot_gaussian_ellipsoid(x(1:2),S(1:2,1:2),sd,p);    
        
        p.axh = axs_v; p.display = 1; 
        plot_gaussian_ellipsoid(x(3:4),S(3:4,3:4),sd,p);
    end
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
    P = P./max(P); 
    
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
function [pt, st] = parse_timing_txt(FILE_PATH)
    fileID = fopen(FILE_PATH, 'r');
    pt = []; st = [];
    
    while ~feof(fileID)
        line = split(fgetl(fileID)); 
        pt(end+1) = str2double(line{3});
        st(end+1) = str2double(line{7});
    end
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%