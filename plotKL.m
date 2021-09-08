% Plot data from multiple diseases output from batchKL.m
clearvars; clc; close all; tic;

% Assumptions and notes
% - assumes datasets for various diseases in Data/[disease name]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       eraged and other state models handle noise

% Directory of some main code and plotting options
thisDir = cd; cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Possible diseases
epiNam = {'marburg', 'mers', 'measles', 'covid', 'ebola'};
% Possible scenarios
scenNam = {'control', 'square-wave', 'cascade', 'boom-bust', 'filtered',...
    'waves', 'noise valley', 'boom-bust-boom', 'rising'};

%% Extract data and plot

% Choose diseases to load
loadDis = [5 4]; nDis = length(loadDis);
loadScen = [2 6]; nScen = length(loadScen);

% All data of interest
data = cell(nDis, nScen);
for i = 1:nDis
    % Disease of interest
    loadname = ['Data/' epiNam{loadDis(i)}]; cd(loadname);
    % Scenarios
    for j = 1:nScen
        data{i, j} = load(['proc1000_' num2str(loadScen(j)) '_' num2str(loadDis(i)) '.mat']);
    end
    cd(thisDir);
    
    %% For each disease generate a plot of both scenarios of interest
    
    % Specific data to plot for first scenario
    D = data{i, 1};

    % Change-points to plot
    figure('Renderer', 'painters', 'Position', [10 10 800 600]);
    subplot(5, 2, [1 3]);
    title(epiNam{loadDis(i)});
    plot(D.tday, D.Rtrue, '-', 'color', 'k', 'linewidth', 2);
    hold on; xlim(D.xlim1);
    plot(D.tday, ones(size(D.Rtrue)), '--', 'color', grey2, 'linewidth', 2);
    plotCIRaw(D.tday', D.Rmf', D.Rlf', D.Ruf', 'b');
    plotCIRaw(D.tday', D.Rm', D.Rl', D.Ru', 'r');
    h = gca; ylim = h.YLim;
    plot([D.pltup D.pltdwn], reshape(repmat(ylim, [1 D.nup+D.ndown]), [2 D.nup+D.ndown]), ...
        '--', 'color', grey2, 'linewidth', 2);
    ylabel('$\hat{R}_s$', 'FontSize', 20);
    grid off; box off; hold off;
    xlim(D.xlim1);
    
    subplot(5, 2, 5);
    plot(D.tday, D.klm, '-', 'color', 'k', 'linewidth', 2);
    hold on; xlim(D.xlim1); h = gca; ylim = h.YLim;
    plot([D.pltup D.pltdwn], reshape(repmat(ylim, [1 D.nup+D.ndown]), [2 D.nup+D.ndown]), ...
        '--', 'color', grey2, 'linewidth', 2);
    grid off; box off; hold off;
    ylabel('D$(p_s | q_s)$', 'FontSize', 20);
    
    subplot(5, 2, [7 9]);
    plot(D.tday, D.Fqm, '-', 'color', 'r', 'linewidth', 2);
    hold on;
    plot(D.tday, D.Fpm, '-', 'color', 'b', 'linewidth', 2);
    plot(D.tday, 0.5*ones(size(D.Rtrue)), '--', 'color', grey2, 'linewidth', 2);
    xlim(D.xlim1); ylim = [0 1.1];
    plot([D.pltup D.pltdwn], reshape(repmat(ylim, [1 D.nup+D.ndown]), [2 D.nup+D.ndown]), ...
        '--', 'color', grey2, 'linewidth', 2);
    hold off; grid off; box off;
    ylabel('P($R_s > 1$)', 'FontSize', 20);
    xlabel('time, $s$ (days)', 'FontSize', 20);
    
    % Specific data to plot for second scenario
    D = data{i, 2};

    % Change-points to plot
    subplot(5, 2, [2 4]);
    plot(D.tday, D.Rtrue, '-', 'color', 'k', 'linewidth', 2);
    hold on; xlim(D.xlim1);
    plot(D.tday, ones(size(D.Rtrue)), '--', 'color', grey2, 'linewidth', 2);
    plotCIRaw(D.tday', D.Rmf', D.Rlf', D.Ruf', 'b');
    plotCIRaw(D.tday', D.Rm', D.Rl', D.Ru', 'r');
    h = gca; ylim = h.YLim;
    plot([D.pltup D.pltdwn], reshape(repmat(ylim, [1 D.nup+D.ndown]), [2 D.nup+D.ndown]), ...
        '--', 'color', grey2, 'linewidth', 2);
    ylabel('$\hat{R}_s$', 'FontSize', 20);
    grid off; box off; hold off;
    xlim(D.xlim1);
    
    subplot(5, 2, 6);
    plot(D.tday, D.klm, '-', 'color', 'k', 'linewidth', 2);
    hold on; xlim(D.xlim1); h = gca; ylim = h.YLim;
    plot([D.pltup D.pltdwn], reshape(repmat(ylim, [1 D.nup+D.ndown]), [2 D.nup+D.ndown]), ...
        '--', 'color', grey2, 'linewidth', 2);
    grid off; box off; hold off;
    ylabel('D$(p_s | q_s)$', 'FontSize', 20);
    
    subplot(5, 2, [8 10]);
    plot(D.tday, D.Fqm, '-', 'color', 'r', 'linewidth', 2);
    hold on;
    plot(D.tday, D.Fpm, '-', 'color', 'b', 'linewidth', 2);
    plot(D.tday, 0.5*ones(size(D.Rtrue)), '--', 'color', grey2, 'linewidth', 2);
    xlim(D.xlim1); ylim = [0 1.1];
    plot([D.pltup D.pltdwn], reshape(repmat(ylim, [1 D.nup+D.ndown]), [2 D.nup+D.ndown]), ...
        '--', 'color', grey2, 'linewidth', 2);
    hold off; grid off; box off;
    ylabel('P($R_s > 1$)', 'FontSize', 20);
    xlabel('time, $s$ (days)', 'FontSize', 20);
    

end