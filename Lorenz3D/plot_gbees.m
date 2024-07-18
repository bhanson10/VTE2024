close all; clc; clear all; 

% plot_gbees.m
% Benjamin Hanson, 2024

%% Define colors
load("colors.mat"); 

%% Initial Condition
const.d = 3; const.T = 1; const.dx=0.5; const.sigma=4; const.b=1; const.r=48;
rv.start=[-11.5; -10; 9.5]; rv.unc = [1; 1; 1];

%% Initializing Figures
initialize_figures(); 

%% Truth
Y0 = rv.start; tspan = [0 50]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,const), tspan, Y0, options);

figure(1);
plot3(Y(:,1),Y(:,2),Y(:,3),'g-','linewidth',1.5,'HandleVisibility','off');  
drawnow;

Y0 = rv.start; tspan = [0 const.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,const), tspan, Y0, options);
plot3(Y(:,1),Y(:,2),Y(:,3),'k-','linewidth',2,'DisplayName','Nominal'); drawnow;

figure(2);
plot3(Y(:,1),Y(:,2),Y(:,3),'k-','linewidth',2,'DisplayName','Nominal'); drawnow;

%{
% Plotting MC trajectories
for k=0:nm
    dir_path = append("./MC/Trajectories/M", num2str(k));
    file_list = dir(fullfile(dir_path, '*.txt'));  % List only .txt files
    num_files = numel(file_list);
    
    x_mc = zeros(size,num_files);
    y_mc = zeros(size,num_files);
    z_mc = zeros(size,num_files);
    for i=0:num_files-1
        file_path = dir_path + "/MC_" + num2str(i) + ".txt"; 
        fileID = fopen(file_path, 'r'); 
        t = str2double(fgetl(fileID)); 
        
        count = 1; 
        while (count < size)
            line = split(fgetl(fileID)); % Read a line as a string
            x_mc(count, i+1) = str2double(line{1});
            y_mc(count, i+1) = str2double(line{2});
            z_mc(count, i+1) = str2double(line{3});
            count = count + 1; 
        end

        % Close the file
        fclose(fileID);
    end

    for i=1:size
        figure(2);
        plot3(x_mc(i,:),y_mc(i,:),z_mc(i,:),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off');
    end
end

Y0 = rv.start; tspan = [0 const.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,const), tspan, Y0, options);

plot3(Y(:,1),Y(:,2),Y(:,3),'k-','linewidth',2,'HandleVisibility','off'); drawnow;

% Plotting MC dots
for k=0:nm
    dir_path = append("./MC/Epochs/M", num2str(k));
    file_list = dir(fullfile(dir_path, '*.txt'));  % List only .txt files
    num_files = numel(file_list); 

    for i=0:num_files-1
        file_path = dir_path + "/MC_" + num2str(i) + ".txt"; 
        fileID = fopen(file_path, 'r'); 
        t = str2double(fgetl(fileID)); 
        
        x_mc = zeros(size); y_mc = zeros(size); z_mc = zeros(size);
        count = 1; 
        while (count < size)
            line = split(fgetl(fileID)); % Read a line as a string
            x_mc(count) = str2double(line{1});
            y_mc(count) = str2double(line{2});
            z_mc(count) = str2double(line{3});
            count = count + 1; 
        end

        figure(2);
        if((k==0)&&(i==0))
            scatter3(x_mc(1),y_mc(1),z_mc(1),2,'k','filled','DisplayName','MC'); 
        end
        scatter3(x_mc,y_mc,z_mc,2,'k','filled','HandleVisibility','off');
        % Close the file
        fclose(fileID);
    end
end
%}

%% GBEES
NM = 1; 
clear p; p.name = "GBEES"; p.alpha = [0.2, 0.4, 0.6]; 
dir_path = "./gbees/v0/Data/PDFs/P";

count = 1;
for nm=0:NM-1

    sub_dir_path = dir_path + num2str(nm); 
    file_list = dir(fullfile(sub_dir_path, '*.txt'));  % List only .txt files
    num_files = numel(file_list);
    
    for i=0:num_files-1
        file_path = sub_dir_path + "/pdf_" + num2str(i) + ".txt";

        [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(file_path);
    
        xest_gbees{count} = zeros(size(x_gbees(1,:)));
        for j=1:n_gbees
            xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
        end

        figure(1); 
        plot_nongaussian_surface(x_gbees,P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
        figure(2); 
        plot_nongaussian_surface(x_gbees,P_gbees,normpdf(3)/normpdf(0), p);
        drawnow; 
        
        count = count + 1;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              FUNCTIONS                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=Lorenz3D(y,const)                          
    f=[const.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -const.b*y(3)+y(1)*y(2)-const.b*const.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures()

    f1 = figure(1); clf; hold all; f1.Position = [150 200 600 475];
    view(-109,14); lighting phong; light('Position',[-1 0 0]); 
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel("x", 'FontSize', 18, 'FontName', 'Times', 'Position',[-10 44 -26]);
    ylabel("y", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 -15 -42]);
    zlabel("z", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 47 8]);
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-20 20])
    xticks([-20 -10 0 10 20])
    xticklabels({'-20','-10','0','10','20'})
    ylim([-30 30])
    yticks([-30 -20 -10 0 10 20 30])
    yticklabels({'-30','-20','-10','0','10', '20', '30'})
    zlim([-30 30])
    zticks([-30 -20 -10 0 10 20 30])
    zticklabels({'-30','-20','-10','0','10', '20', '30'})
    
    f2 = figure(2); clf; hold all; f2.Position = [750 200 600 475];
    view(-109,14); lighting phong; light('Position',[-1 0 0]);
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel("x", 'FontSize', 18, 'FontName', 'Times', 'Position',[-10 44 -26]);
    ylabel("y", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 -15 -42]);
    zlabel("z", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 47 8]);
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-20 20])
    xticks([-20 -10 0 10 20])
    xticklabels({'-20','-10','0','10','20'})
    ylim([-30 30])
    yticks([-30 -20 -10 0 10 20 30])
    yticklabels({'-30','-20','-10','0','10', '20', '30'})
    zlim([-30 30])
    zticks([-30 -20 -10 0 10 20 30])
    zticklabels({'-30','-20','-10','0','10', '20', '30'})
    set(gca, 'FontName' , 'Times');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        x(count, :) = [str2double(line{2});str2double(line{3});str2double(line{4})];
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%