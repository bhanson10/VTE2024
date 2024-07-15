clear all; close all; clc; 

DATASET = "./Data"; TU = 1;
name1 = 'C w/ Python extension';
name2 = 'C++';
time_unit = '(TU)';

f1 = fopen(DATASET + "/runtime_extension.txt", 'r');
f2 = fopen(DATASET + "/runtime_c.txt", 'r');


pf = 1; % Choose which speed test to perform, 1 or 2

if(pf==1)
    [rt1, rt2, pt1, pt2, s1, s2] = parse_files_1(f1,f2);

    f1 = figure(1); clf; f1.Position = [200 350 700 475]; ax = axes; 
    l = legend; l.Location = "Northwest"; l.FontSize = 14; l.FontName = "Times"; 
    set(ax, 'FontName' , 'Times','FontSize',14);
    xlabel("Simulation time " + time_unit, 'FontSize', 18, 'FontName', 'Times');
    xlim([rt1(1).*TU rt1(end).*TU]);

    yyaxis left; hold on; 
    ax.YAxis(1).Color = [0 0 1];
    plot(rt1.*TU, pt1, 'b-', 'LineWidth', 1, 'DisplayName',name1);
    plot(rt2.*TU, pt2, 'b--', 'LineWidth', 1,'DisplayName',name2);
    ylabel("Program time (s)", 'FontSize', 18, 'FontName', 'Times');
    
    yyaxis right; hold on; 
    ax.YAxis(2).Color = [1 0 0];
    plot(rt1.*TU, s1, 'r-', 'LineWidth', 1, 'DisplayName',name1);
    plot(rt2.*TU, s2, 'r--', 'LineWidth', 1, 'DisplayName',name2);
    ylabel("Number of cells", 'FontSize', 18, 'FontName', 'Times');
    
    % Normalized time
    f2 = figure(2); clf; f2.Position = [900 350 700 475]; ax = axes; 
    l = legend; l.Location = "Northeast"; l.FontSize = 14; l.FontName = "Times"; 
    set(ax, 'FontName' , 'Times','FontSize',14);
    xlabel("Simulation time " + time_unit, 'FontSize', 18, 'FontName', 'Times');
    xlim([rt1(2).*TU rt1(end).*TU]);

    yyaxis left; hold on; 
    pt_norm = pt1./pt2;
    plot(rt1(2:end).*TU, pt_norm(2:end), 'b-', 'LineWidth', 1, 'DisplayName', append(name1, ' / ', name2));
    ylabel("Normalized program time", 'FontSize', 18, 'FontName', 'Times');
    
    yyaxis right; hold on; 
    ylim([0.5,1.5])
    s_norm = s1./s2;
    plot(rt1(2:end).*TU, s_norm(2:end), 'r-', 'LineWidth', 1, 'DisplayName', append(name1, ' / ', name2));
    ylabel("Normalized cell number", 'FontSize', 18, 'FontName', 'Times');

    ax.YAxis(1).Color = [0 0 1];
    ax.YAxis(2).Color = [1 0 0];

elseif(pf==2)
    [rt1, rt2, pt1, pt2, s1] = parse_files_2(f1,f2);

    f1 = figure(1); clf; hold all; f1.Position = [400 250 700 475]; 
    l = legend; l.Location = "Northwest"; l.FontSize = 14; l.FontName = "Times"; 
    set(gca, 'FontName' , 'Times','FontSize',14);
    xlabel("Simulation time " + time_unit, 'FontSize', 18, 'FontName', 'Times');
    ylabel("Program time (s)", 'FontSize', 18, 'FontName', 'Times');
    xlim([rt1(1) rt1(end)]);

    yyaxis left;
    plot(rt1, pt1, 'b-', 'LineWidth', 1, 'DisplayName',name1);
    plot(rt2, pt2, 'b--', 'LineWidth', 1,'DisplayName',name2);

    yyaxis right;
    plot(rt1, s1, 'r-', 'LineWidth', 1, 'DisplayName',name1);
    ylabel("Number of cells", 'FontSize', 18, 'FontName', 'Times');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rt1, rt2, pt1, pt2, s1, s2] = parse_files_1(f1,f2)
    rt1 = []; rt2 = rt1; pt1 = rt1; pt2 = rt1; s1 = rt1; s2 = rt1; 
    
    while ~feof(f1)
        line = split(fgetl(f1)); 
        pt1(end+1) = str2double(line{5});
        rt1(end+1) = str2double(line{9});
        size = strsplit(line{13},'/');
        s1(end+1) = str2double(size{1}); 
    end

    while ~feof(f2)
        line = split(fgetl(f2)); 
        pt2(end+1) = str2double(line{5});
        rt2(end+1) = str2double(line{9});
        size = strsplit(line{13},'/');
        s2(end+1) = str2double(size{1}); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rt1, rt2, pt1, pt2, s1] = parse_files_2(f1,f2)
    rt1 = []; rt2 = rt1; pt1 = rt1; pt2 = rt1; s1 = rt1; 
    
    while ~feof(f1)
        line = split(fgetl(f1)); 
        pt1(end+1) = str2double(line{5});
        rt1(end+1) = str2double(line{9});
        size = strsplit(line{13},'/');
        s1(end+1) = str2double(size{1}); 
    end

    while ~feof(f2)
        line = split(fgetl(f2)); 
        pt2(end+1) = str2double(line{5});
        rt2(end+1) = str2double(line{9});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%