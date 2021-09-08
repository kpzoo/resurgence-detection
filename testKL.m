% Test real-time behaviour with KL statistics
clearvars; clc; close all; tic;

% Assumptions and notes
% - test real time via filtering, smoothing and sequential smoothing
% - examine weekly averaged and other state models handle noise

% Directory and where saving
thisDir = cd; saveFol = 'Results/examples'; 
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Simulate epidemic 

% Define horizon for prediction
horiz = 30; disp(['Predicting ' num2str(horiz) ' days ahead']);
% Define number of days to simulate
nday0 = 301; tday0 = 1:nday0;
% Choose scenarios and fix serial interval
scenNo = 9; epiNo = 4;

% Save names
switch(epiNo)
    case 1
        saveFol = 'Data/marburg'; 
    case 2
        saveFol = 'Data/mers';
    case 3
        saveFol = 'Data/measles';
    case 4
        saveFol = 'Data/covid'; 
    case 5
        saveFol = 'Data/ebola'; 
end

% Simulation parameters and warning
Iwarn = 1; simVals = setupScenario(scenNo); simVals.offset = 0;
% Simulate epidemic scenarios and truncate
while Iwarn
    [Iday, Lday, Rtrue, tday, Iwarn, distvals] = epiSimScenRealTime(scenNo, nday0, epiNo, simVals);
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Truncated observation period 
nday = length(tday);

% Saving data and figs
namstr = [num2str(nday) '_' num2str(horiz) '_' num2str(epiNo) '_'];

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

%% Estimate Rt from standard state model 

% Standard EpiFilter and EpiSmoother estimates
[~, Rlowf, Rhighf, Rmeanf, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Lday, Iday);
[~, Rlow, Rhigh, Rmean, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% Standard EpiSmoother one-step-ahead predictions 
[Imean, predInt] = recursPredict(Rgrid, qR, Lday, Rmean, max(Iday));
Ilow = predInt(:,1)'; Ihigh = predInt(:,2)';
[Imeanf, predIntf] = recursPredict(Rgrid, pR, Lday, Rmean, max(Iday));
Ilowf = predIntf(:,1)'; Ihighf = predIntf(:,2)';

% Main horizon function for statistics on Rt and It (standard)
[Rhoriz, Ihoriz, Rstats, Istats] = predHoriz(Rgrid, m, eta, horiz, qR(end, :), Lday, Iday, p0, distvals, nday);
[Rhorizf, Ihorizf, Rstatsf, Istatsf] = predHoriz(Rgrid, m, eta, horiz, pR(end, :), Lday, Iday, p0, distvals, nday);

% Twice the largest width of fitted R
Rw_max = max(Rhigh - Rlow)*ones(1, horiz);
% Width of R from horizon 
Rw_horiz = Rhoriz(:,3) - Rhoriz(:,1);
Rw_horizf = Rhorizf(:,3) - Rhorizf(:,1);

% Compute Kullback Liebler divergence
kl = zeros(1, nday); 
for i = 1:nday
   pi = pR(i, :); qi = qR(i, :);
   % Elements of KL vector
   klvec = pi.*log(pi./qi);
   % Remove edge cases
   klvec(isinf(klvec) | isnan(klvec)) = 0;
   kl(i) = sum(klvec);
end

% Compute Jensen Shannon divergence
js = zeros(1, nday); 
for i = 1:nday
   pi = pR(i, :); qi = qR(i, :);
   mi = 0.5*(pi + qi);
   % Elements of JS vector
   jsvec = pi.*log(pi./mi) + qi.*log(qi./mi);
   % Remove edge cases
   jsvec(isinf(jsvec) | isnan(jsvec)) = 0;
   js(i) = 0.5*sum(jsvec);
end

% CDFs of R > 1
id1 = find(Rgrid <= 1, 1, 'last');
Fq = 1 - sum(qR(:, 1:id1), 2)';
Fp = 1 - sum(pR(:, 1:id1), 2)';

%% Visualisation of results

% Axes handles and horizon times
ax1 = zeros(1, 2); ax2 = ax1;
thoriz = tday(end)+1:tday(end)+horiz;

% Divergences against true R
figure;
%plot(tday, kl, '-', 'color', 'r', 'linewidth', 2);
hold on;
plot(tday, js, '-', 'color', 'b', 'linewidth', 2);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
xlim([tday(10), tday(end)]);
ylabel('JS$(p | q)$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);


% CDFs of R > 1 against true R
figure;
plot(tday, Fq, '-', 'color', 'r', 'linewidth', 2);
hold on;
plot(tday, Fp, '-', 'color', 'b', 'linewidth', 2);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
xlim([tday(10), tday(end)]);
ylabel('P$(R > 1)$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);

% Comparison R estimates and I predictions on known data
figure;
ax1(1) = subplot(2, 1, 1);
plot(tday(2:end), Rtrue(2:end), '-', 'color', 'k', 'linewidth', 2);
hold on; xlim([tday(2) tday(end)]);
plotCIRaw(tday', Rmean', Rlow', Rhigh', 'r');
grid off; box off; hold off;
ylabel('$\hat{R}_{s}$', 'FontSize', 20);
ax1(2) = subplot(2, 1, 2);
plot(tday(2:end)', Iday(2:end)', '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim([tday(2) tday(end)]);
plotCIRaw(tday(2:end)', Imean', Ilow', Ihigh', 'b');
grid off; box off; hold off;
ylabel('$\hat{I}_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
linkaxes(ax1, 'x'); 


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
plot(tday(2:end)', Iday(2:end)', '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim([tday(2) thoriz(end)]);
plotCIRaw(tday(2:end)', Imean', Ilow', Ihigh', 'b');
plotCIRaw(thoriz', Istats.mean, Ihoriz(:, 1), Ihoriz(:, 3), 'r'); 
grid off; box off; hold off;
ylabel('$\hat{I}_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
linkaxes(ax1, 'x'); 

% Comparison R estimates and I predictions on known data (modified)
figure;
ax1(1) = subplot(2, 1, 1);
plot(tday(2:end), Rtrue(2:end), '-', 'color', 'k', 'linewidth', 2);
hold on; xlim([tday(2) tday(end)]);
plotCIRaw(tday', Rmeanf', Rlowf', Rhighf', 'r');
grid off; box off; hold off;
ylabel('$\hat{R}_{s}$', 'FontSize', 20);
ax1(2) = subplot(2, 1, 2);
plot(tday(2:end)', Iday(2:end)', '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim([tday(2) tday(end)]);
plotCIRaw(tday(2:end)', Imeanf', Ilowf', Ihighf', 'b');
grid off; box off; hold off;
ylabel('$\hat{I}_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
linkaxes(ax1, 'x'); 


% Include predictions of R and I across horizon (modified)
figure;
ax1(1) = subplot(2, 1, 1);
plot(tday(2:end), Rtrue(2:end), '-', 'color', 'k', 'linewidth', 2);
hold on; xlim([tday(2) thoriz(end)]);
plotCIRaw(tday', Rmeanf', Rlowf', Rhighf', 'r');
plotCIRaw(thoriz', Rstatsf.mean, Rhorizf(:, 1), Rhorizf(:, 3), 'r'); 
grid off; box off; hold off;
ylabel('$\hat{R}_{s}$', 'FontSize', 20);
ax1(2) = subplot(2, 1, 2);
plot(tday(2:end)', Iday(2:end)', '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim([tday(2) thoriz(end)]);
plotCIRaw(tday(2:end)', Imeanf', Ilowf', Ihighf', 'b');
plotCIRaw(thoriz', Istatsf.mean, Ihorizf(:, 1), Ihorizf(:, 3), 'r'); 
grid off; box off; hold off;
ylabel('$\hat{I}_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
linkaxes(ax1, 'x'); 

% Direct comparison of R estimates from pR and qR
figure;
plot(tday(2:end), Rtrue(2:end), '-', 'color', 'k', 'linewidth', 2);
hold on; xlim([tday(2) tday(end)]);
plotCIRaw(tday', Rmean', Rlow', Rhigh', 'r');
plotCIRaw(tday', Rmeanf', Rlowf', Rhighf', 'b');
grid off; box off; hold off;
ylabel('$\hat{R}_{s}$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);