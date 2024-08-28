clear all; close all; clc;

% synoydic_to_inertial.m
% Benjamin Hanson, 2024

%% colors
C = [238, 102, 119;  % Red
     68,  119, 170;  % Blue
     255, 140, 0;    % Orange
     34,  139, 34;]; % Green

C  = C/255;

%% user input
sys = "SaEn";   % three-body system
orb = "DPO";    % family of periodic orbits
id = 1176;      % orbit identifier

%% initialize system properties
ic = read_ic(sys, orb, id);
[prop.mu, prop.LU, prop.TU, prop.pri_r, prop.sec_r] = get_prop(sys);
prop.U = [prop.LU, prop.LU, prop.LU/prop.TU, prop.LU/prop.TU]; 
models = {'mc', 'ukf', 'bpf', 'gbees'};
initialize_figures(prop);

%% nominal trajectory
rv_DPO_synodic = parse_monte_txt("rv_dpo_synodic.txt");
rv_DPO_inertial = parse_monte_txt("rv_dpo_inertial.txt");
rv_enc_inertial = parse_monte_txt("rv_enc_inertial.txt");
rv_DPO_synodic = rv_DPO_synodic.*prop.U;

figure(1); 
nexttile(1); 
plot(rv_DPO_synodic(:,1),rv_DPO_synodic(:,2),'k-','LineWidth', 2,'HandleVisibility','off');
drawnow;
nexttile(2); 
plot(rv_DPO_synodic(:,3),rv_DPO_synodic(:,4),'k-','LineWidth', 2,'HandleVisibility','off');
drawnow;

clear L; clear LH; 
LH(1) = plot(NaN,NaN,'o','MarkerSize', 10, 'MarkerFaceColor',[1,167/255,254/255],'MarkerEdgeColor','none');
L{1} = "Enceladus\,\,\,";
LH(2) = plot(NaN,NaN,'k-', 'LineWidth',2);
L{2} = "DPO\,\,\,";
leg = legend(LH, L, 'Orientation', 'Horizontal', 'FontSize', 18, 'FontName', 'times', 'Interpreter', 'latex');
leg.Layout.Tile = 'south';
drawnow; 

figure(2); 
nexttile(1); 
plot(rv_DPO_inertial(:,1),rv_DPO_inertial(:,2),'k-','LineWidth', 2,'HandleVisibility','off');
plot(rv_enc_inertial(:,1),rv_enc_inertial(:,2),'m-','LineWidth', 2,'HandleVisibility','off');
drawnow;
nexttile(2); 
plot(rv_DPO_inertial(:,3),rv_DPO_inertial(:,4),'k-','LineWidth', 2,'HandleVisibility','off');
plot(rv_enc_inertial(:,3),rv_enc_inertial(:,4),'m-','LineWidth', 2,'HandleVisibility','off');
drawnow;

clear L; clear LH; 
LH(1) = plot(NaN,NaN,'o','MarkerSize', 10, 'MarkerFaceColor',[232/255,230/255,200/255],'MarkerEdgeColor','none');
L{1} = "Saturn\,\,\,";
LH(2) = plot(NaN,NaN,'m-', 'LineWidth',2);
L{2} = "Enceladus\,\,\,";
LH(3) = plot(NaN,NaN,'k-', 'LineWidth',2);
L{3} = "DPO\,\,\,";
leg = legend(LH, L, 'Orientation', 'Horizontal', 'FontSize', 18, 'FontName', 'times', 'Interpreter', 'latex');
leg.Layout.Tile = 'south';
drawnow; 

figure(3); 
nexttile(1); 
plot(rv_DPO_synodic(:,1),rv_DPO_synodic(:,2),'k-','LineWidth', 2,'HandleVisibility','off');
drawnow;
nexttile(2); 
plot(rv_DPO_inertial(:,1),rv_DPO_inertial(:,2),'k-','LineWidth', 2,'HandleVisibility','off');
plot(rv_enc_inertial(:,1),rv_enc_inertial(:,2),'m-','LineWidth', 2,'HandleVisibility','off');
drawnow;

clear L; clear LH; 
LH(1) = plot(NaN,NaN,'o','MarkerSize', 10, 'MarkerFaceColor',[1,167/255,254/255],'MarkerEdgeColor','none');
L{1} = "Enceladus\,\,\,";
LH(2) = plot(NaN,NaN,'o','MarkerSize', 10, 'MarkerFaceColor',[232/255,230/255,200/255],'MarkerEdgeColor','none');
L{2} = "Saturn\,\,\,";
LH(3) = plot(NaN,NaN,'k-', 'LineWidth',2);
L{3} = "DPO\,\,\,";
leg = legend(LH, L, 'Orientation', 'Horizontal', 'FontSize', 18, 'FontName', 'times', 'Interpreter', 'latex');
leg.Layout.Tile = 'south';
drawnow; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ic = read_ic(sys, orb, id)
    file_str = "./ephem/" + sys + "_" + orb + "_" + num2str(id) + ".csv";
    ic_table = readmatrix(file_str);
    ic_data = ic_table(1,:); 
    ic.id = ic_data(1); 
    ic.state = [ic_data(2); ic_data(3); ic_data(5); ic_data(6)];
    ic.J = ic_data(8);
    ic.T = ic_data(9); 
    ic.SI = ic_data(11);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu, LU, TU, pri_r, sec_r] = get_prop(sys)
    switch sys
        case 'SuEa'
            mu = 3.054200000000000E-6;
            LU = 149597871;
            TU = 5022635;
            pri_r = 696340;
            sec_r = 6378.0;
        case 'SuMa'
            mu = 3.227154996101724E-7;
            LU = 208321282;
            TU = 8253622;
            pri_r = 696340;
            sec_r = 3389.5;
        case 'EaMo'
            mu = 1.215058560962404E-2;
            LU = 389703;
            TU = 382981;
            pri_r = 6378.14;
            sec_r = 1737.1;
        case 'MaPh'
            mu = 1.611081404409632E-8;
            LU = 9468;
            TU = 4452;
            pri_r = 3389.5;
            sec_r = 11.3;
        case 'JuEu'
            mu = 2.528017528540000E-5;
            LU = 668519;
            TU = 48562;
            pri_r = 69911;
            sec_r = 1560.8;
        case 'SaEn'
            mu = 1.901109735892602E-7;
            LU = 238529;
            TU = 18913;
            pri_r = 58232;
            sec_r = 252.1;
        case 'SaTi'
            mu = 2.366393158331484E-4;
            LU = 1195677;
            TU = 212238;
            pri_r = 58232;
            sec_r = 2574.7;
        otherwise
            error("System is not in database. Make sure case is correct.")
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function axs = initialize_figures(prop)
    f1 = figure(1); clf; hold all; f1.Position = [200 100 1200 700];
    tiledlayout(1, 2, 'TileSpacing','compact');
    
    axs{1} = nexttile(1); hold all;
    axis("equal")
    set(gca, 'FontName', 'Times', 'FontSize', 14);
    xlabel("synodic $x$ (km)", "Interpreter","latex")
    ylabel("synodic $y$ (km)", "Interpreter","latex")
    annotation('arrow', [0.413, 0.413 - 0.03], [0.2, 0.2], 'Color', 'black', 'LineWidth', 1);
    text(2.3903E5, -980, 'to Saturn', 'FontName', 'times', 'Color', 'black', 'FontSize', 12);
    p2 = nsidedpoly(1000, 'Center', [(1-prop.mu)*prop.LU, 0], 'Radius', prop.sec_r);
    plot(p2, 'FaceColor', 'm', 'EdgeColor','none');
    
    nexttile(2); hold all;
    axis("equal")
    set(gca, 'FontName', 'Times', 'FontSize', 14);
    xlabel("synodic $\dot{x}$ (km/s)", "Interpreter","latex")
    ylabel("synodic $\dot{y}$ (km/s)", "Interpreter","latex")

    f2 = figure(2); clf; hold all; f2.Position = [200 100 1200 700];
    tiledlayout(1, 2, 'TileSpacing','compact');
    
    axs{2} = nexttile(1); hold all;
    axis("equal")
    set(gca, 'FontName', 'Times', 'FontSize', 14);
    xlabel("inertial $X$ (km)", "Interpreter","latex")
    ylabel("inertial $Y$ (km)", "Interpreter","latex")
    p1 = nsidedpoly(1000, 'Center', [0, 0], 'Radius', prop.pri_r);
    plot(p1, 'FaceColor', [189/255, 183/255, 107/255], 'EdgeColor','none');
    
    nexttile(2); hold all;
    axis("equal")
    set(gca, 'FontName', 'Times', 'FontSize', 14);
    xlabel("inertial $\dot{X}$ (km/s)", "Interpreter","latex")
    ylabel("inertial $\dot{Y}$ (km/s)", "Interpreter","latex")

    f3 = figure(3); clf; hold all; f3.Position = [200 100 1200 700];
    tiledlayout(1, 2, 'TileSpacing','compact');

    axs{3} = nexttile(); hold all;
    axis("equal")
    set(gca, 'FontName', 'Times', 'FontSize', 14);
    xlabel("synodic $x$ (km)", "Interpreter","latex")
    ylabel("synodic $y$ (km)", "Interpreter","latex")
    annotation('arrow', [0.414, 0.414 - 0.03], [0.2, 0.2], 'Color', 'black', 'LineWidth', 1);
    text(2.3903E5, -980, 'to Saturn', 'FontName', 'times', 'Color', 'black', 'FontSize', 12);
    p2 = nsidedpoly(1000, 'Center', [(1-prop.mu)*prop.LU, 0], 'Radius', prop.sec_r);
    plot(p2, 'FaceColor', 'm', 'EdgeColor','none');

    nexttile(2); hold all;
    axis("equal")
    set(gca, 'FontName', 'Times', 'FontSize', 14);
    xlabel("inertial $X$ (km)", "Interpreter","latex")
    ylabel("inertial $Y$ (km)", "Interpreter","latex")
    p1 = nsidedpoly(1000, 'Center', [0, 0], 'Radius', prop.pri_r);
    plot(p1, 'FaceColor',  [189/255, 183/255, 107/255], 'EdgeColor','none');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = PCR3BP(t, x, prop)
    x1 = [x(3); x(4); 2*x(4)+x(1)-(prop.mu*(x(1)-1+prop.mu)/(((x(1)-1+prop.mu)^2+x(2)^2)^(1.5)))-((1-prop.mu)*(x(1)+prop.mu)/(((x(1)+prop.mu)^2+x(2)^2)^(1.5))); -2*x(3)+x(2)-(prop.mu*x(2)/(((x(1)-1+prop.mu)^2+x(2)^2)^(1.5)))-((1-prop.mu)*x(2)/(((x(1)+prop.mu)^2+x(2)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv = parse_monte_txt(filename)
    fileID = fopen(filename, 'r');
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        rv(count, :) = [str2double(line{1}); str2double(line{2}); str2double(line{4}); str2double(line{5})];
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%