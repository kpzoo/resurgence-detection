% Test prediction horizon of smoothed distributions
clearvars; clc; close all; tic;

% Assumptions and notes
% - examine how the posterior over R widens with time and no data
% - recursively transform this into predictions of incidence

% Directory and where saving
thisDir = cd; saveFol = 'Results/'; 
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Simulate epidemic and estimate Rt

% Define horizon for prediction
horiz = 30; disp(['Predicting ' num2str(horiz) ' days ahead']);
% Define number of days to simulate
nday0 = 301; tday0 = 1:nday0;
% Choose scenarios and fix serial interval
scenNo = 7; distNo = 2;

% Simulation parameters and warning
Iwarn = 1; simVals = setupScenario(scenNo); simVals.offset = 0;
% Simulate epidemic scenarios and truncate
while Iwarn
    [Iday, Lday, Rtrue, tday, Iwarn, distvals] = ...
        epiSimScenOffset(scenNo, nday0, distNo, simVals);
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Truncated observation period - same for all demes
nday = length(tday);

% Saving data and figs
namstr = [num2str(nday) '_' num2str(horiz) '_' num2str(distNo) '_'];

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% EpiFilter and EpiSmoother estimates
[~, ~, ~, ~, pR, pRup, pstate] = runEpiFilterState(Rgrid, m, eta, nday, p0, Lday, Iday);
[~, Rlow, Rhigh, Rmean, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiSmoother one-step-ahead predictions 
[Imean, predInt] = recursPredict(Rgrid, qR, Lday, Rmean, 2*max(Iday));
Ilow = predInt(:,1)'; Ihigh = predInt(:,2)';

%% Prediction of Rt and It

% Main horizon function for statistics on Rt and It
[Rhoriz, Ihoriz, Rstats, Istats] = predHoriz(Rgrid, m, eta, horiz, qR(end, :), Lday, Iday, p0, distvals, nday);

% Twice the largest width of fitted R
Rw_max = max(Rhigh - Rlow)*ones(1, horiz);
% Width of R from horizon 
Rw_horiz = Rhoriz(:,3) - Rhoriz(:,1);
% Time of crossover
tcross = find(Rw_horiz >= Rw_max(1), 1, 'first');
disp(['Max horizon is ' num2str(tcross)]);

%% Visualisation of results

% Axes handles and horizon times
ax1 = zeros(1, 2); ax2 = ax1;
thoriz = tday(end)+1:tday(end)+horiz;

% Comparison R estimates and I predictions on known data
figure;
ax1(1) = subplot(2, 1, 1);
plot(tday(2:end), Rtrue(2:end), '-', 'color', 'k', 'linewidth', 2);
hold on; xlim([tday(2) tday(end)]);
plotCIRaw(tday', Rmean', Rlow', Rhigh', 'r');
grid off; box off; hold off;
ylabel('$\hat{R}_{s}$', 'FontSize', 20);
ax1(2) = subplot(2, 1, 2);
% scatter(tday(2:end)', Iday(2:end)', 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
%     'k', 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
plot(tday(2:end)', Iday(2:end)', '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim([tday(2) tday(end)]);
plotCIRaw(tday(2:end)', Imean', Ilow', Ihigh', 'b');
grid off; box off; hold off;
ylabel('$\hat{I}_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
linkaxes(ax1, 'x'); 
if saveFig
    cd(saveFol);
    saveas(gcf, ['est_' namstr 'known'], 'fig');
    cd(thisDir);
end

% Include predictions of R and I across horizon
figure;
ax1(1) = subplot(2, 1, 1);
plot(tday(2:end), Rtrue(2:end), '-', 'color', 'k', 'linewidth', 2);
hold on; xlim([tday(2) thoriz(end)]);
plotCIRaw(tday', Rmean', Rlow', Rhigh', 'r');
plotCIRaw(thoriz', Rstats.mean, Rhoriz(:, 1), Rhoriz(:, 3), 'r'); 
grid off; box off; hold off;
ylabel('$\hat{R}_{s}$', 'FontSize', 20);
ax1(2) = subplot(2, 1, 2);
% scatter(tday(2:end)', Iday(2:end)', 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
%     'k', 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
plot(tday(2:end)', Iday(2:end)', '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim([tday(2) thoriz(end)]);
plotCIRaw(tday(2:end)', Imean', Ilow', Ihigh', 'b');
plotCIRaw(thoriz', Istats.mean, Ihoriz(:, 1), Ihoriz(:, 3), 'r'); 
grid off; box off; hold off;
ylabel('$\hat{I}_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
linkaxes(ax1, 'x'); 
if saveFig
    cd(saveFol);
    saveas(gcf, ['est_' namstr 'horiz'], 'fig');
    cd(thisDir);
end

figure;
plot(thoriz - tday(end), Rw_horiz, 'LineWidth', 2);
hold on;
plot(thoriz - tday(end), Rw_max, 'k--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\Delta w(R)_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
if saveFig
    cd(saveFol);
    saveas(gcf, ['tcross_' namstr], 'fig');
    cd(thisDir);
end