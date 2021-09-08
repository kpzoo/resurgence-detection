% Examine delay in real time behaviour
clearvars; clc; close all; tic;

% Assumptions and notes
% - compare CDFs from filtering, smoothing and sequential smoothing
% - get delay due to lack of knowledge and average over runs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         eraged and other state models handle noise

% Directory and where saving
thisDir = cd; saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Simulate epidemic 

% Set number of epidemics 
M = 1000; disp(['Simulating ' num2str(M) ' epidemics']);

% Define number of days to simulate
nday0 = 301; tday0 = 1:nday0;
% Choose scenarios and fix serial interval
scenNo = 2; epiNo = 5;

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

% Epidemics variables to store
Iday = cell(1, M); Lday = Iday; tday = Iday;
% Simulate M epidemic replicates
for i = 1:M
    % Simulation parameters and warning
    Iwarn = 1; simVals = setupScenario(scenNo); simVals.offset = 0;
    while Iwarn
        [Iday{i}, Lday{i}, Rtrue, tday{i}, Iwarn, distvals] = ...
            epiSimScenRealTime(scenNo, nday0, epiNo, simVals);
        if scenNo == 2 || scenNo == 6
            % Remove sequences with few cases at end
            idrem = sum(Iday{i}(end-5:end)) < 20;
            if idrem
                Iwarn = 1;
            end
        end
    end
    if Iwarn
        warning('Sequences of zero incidence');
    end
end
% Truncated observation period 
nday = cellfun(@length, tday); clc;
if length(unique(nday)) == 1
    % All incidence curves should be same length
    nday = nday(1); tday = tday{1};
else 
    error('Incidence curves inconsistent');
end

% Find crossing points of Rtrue
tup = cell(1, 1); tdown = tup; jup = 0; jdown = 0;
for i = 2:nday
    % Find resurgence crossings
    if sign(Rtrue(i-1) - 1) <= 0 && sign(Rtrue(i) - 1) > 0
        jup = jup + 1; tup{jup} = i;
    end
    % Find crossings indicating control
    if sign(Rtrue(i-1) - 1) >= 0 && sign(Rtrue(i) - 1) < 0
        jdown = jdown + 1; tdown{jdown} = i;
    end
end
tup = cell2mat(tup); tdown = cell2mat(tdown);
nup = length(tup); ndown = length(tdown);

% Saving data and figs
namstr = [num2str(M) '_' num2str(scenNo) '_' num2str(epiNo)];

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

%% Estimate Rt and get CDFs

% R and incidence variables
Rmean = cell(1, M); Rhigh = Rmean; Rlow = Rmean;
Rmeanf = Rmean; Rhighf = Rmean; Rlowf = Rmean;

% KL divergence and CDFs
kl = zeros(nday, M); js = kl; Fq = kl; Fp = kl;

for i = 1:M
    % Standard EpiFilter and EpiSmoother estimates
    [~, Rlowf{i}, Rhighf{i}, Rmeanf{i}, pR, pRup, pstate] = ...
        runEpiFilter(Rgrid, m, eta, nday, p0, Lday{i}, Iday{i});
    [~, Rlow{i}, Rhigh{i}, Rmean{i}, qR] = ...
        runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);
    
    % Compute KL divergence
    for j = 1:nday
        pi = pR(j, :); qi = qR(j, :);
        % Elements of KL vector
        klvec = pi.*log(pi./qi);
        % Remove edge cases
        klvec(isinf(klvec) | isnan(klvec)) = 0;
        kl(j, i) = sum(klvec);
    end
    
    % Compute Jensen Shannon divergence
    js = zeros(1, nday);
    for j = 1:nday
        pi = pR(j, :); qi = qR(j, :);
        mi = 0.5*(pi + qi);
        % Elements of JS vector
        jsvec = pi.*log(pi./mi) + qi.*log(qi./mi);
        % Remove edge cases
        jsvec(isinf(jsvec) | isnan(jsvec)) = 0;
        js(j, i) = 0.5*sum(jsvec);
    end
    
    % Complementary CDFs for R > 1
    id1 = find(Rgrid <= 1, 1, 'last');
    Fq(:, i) = 1 - sum(qR(:, 1:id1), 2)';
    Fp(:, i) = 1 - sum(pR(:, 1:id1), 2)';
    
    disp(['Completed ' num2str(i) ' of ' num2str(M)]);
end

% Mean and quantiles
Fqm = mean(Fq, 2); Fpm = mean(Fp, 2); klm = mean(kl, 2); jsm = mean(js, 2);
Fqq = quantile(Fq, [0.025 0.975], 2)'; Fpq = quantile(Fp, [0.025 0.975], 2)';
% Variances of complementary CDFs
Fqvar = var(Fq, [], 2); Fpvar = var(Fp, [], 2); 

%% Examine timing delays in signalling changes around 1

% Probability countours for checking F(1) and 1 - F(1)
pr = 0.5:0.2:0.9; lenr = length(pr);

% Delays for upward crossings
tqset = zeros(nup, M, lenr); tpset = tqset;
for k = 1:nup
    % Do delay for each run
    for j = 1:M
        for i = 1:lenr
            tqtemp = find(Fq(:, j) >= pr(i)); tqtemp = tqtemp(tqtemp > tup(k)-5);
            tptemp = find(Fp(:, j) >= pr(i)); tptemp = tptemp(tptemp > tup(k)-5);
            tqset(k, j, i) = tqtemp(1); tpset(k, j, i) = tptemp(1);
        end
    end
end

% Different between qR and pR based complementary CDFs
tdiffup = tpset - tqset; tmeanup = zeros(nup, lenr);
for i = 1:nup
    ttemp = squeeze(tdiffup(i, :, :));
    % Delay relative to mean serial interval
    tmeanup(i, :) = mean(ttemp)/distvals.omega;
end

% Delays for downward crossings
tqset2 = zeros(ndown, M, lenr); tpset2 = tqset2;
for k = 1:ndown
    % Do delay for each run
    for j = 1:M
        for i = 1:lenr
            tqtemp = find(1-Fq(:, j) >= pr(i)); tqtemp = tqtemp(tqtemp > tdown(k)-5);
            tptemp = find(1-Fp(:, j) >= pr(i)); tptemp = tptemp(tptemp > tdown(k)-5);
            tqset2(k, j, i) = tqtemp(1); tpset2(k, j, i) = tptemp(1);
        end
    end
end

% Different between qR and pR based complementary CDFs
tdiffdown = tpset2 - tqset2; tmeandown = zeros(ndown, lenr);
for i = 1:ndown
    ttemp = squeeze(tdiffdown(i, :, :));
    % Delay relative to mean serial interval
    tmeandown(i, :) = mean(ttemp)/distvals.omega;
end

% Envelope of mean estimates
Rmean = cell2mat(Rmean'); Rmeanf = cell2mat(Rmeanf');
Rm = mean(Rmean); Rmf = mean(Rmeanf);
Rl = quantile(Rmean, 0.025); Rlf = quantile(Rmeanf, 0.025);
Ru = quantile(Rmean, 1-0.025); Ruf = quantile(Rmeanf, 1-0.025);

%% Visualisation of results

% Limits for axes
xlim1 = [tday(20) tday(end)];

% CDFs with confidence intervals
figure;
yyaxis right
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
ylabel('$R_{s}$', 'FontSize', 20);
yyaxis left
plotCIRaw(tday', Fqm, Fqq(1, :)', Fqq(2, :)', 'r');
hold on;
plotCIRaw(tday', Fpm, Fpq(1, :)', Fpq(2, :)', 'b');
ylabel('$1 - F_s(1)$', 'FontSize', 20);
grid off; box off; hold off;
xlabel('time, $s$ (days)', 'FontSize', 20);

% KL divergence against true R
figure;
yyaxis right
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
ylabel('$R_{s}$', 'FontSize', 20);
yyaxis left
plot(tday, klm, '-', 'color', 'r', 'linewidth', 2);
hold on;
plot(tday, jsm, '-', 'color', 'b', 'linewidth', 2);
ylabel('KL$(p | q)$', 'FontSize', 20);
grid off; box off; hold off;
xlabel('time, $s$ (days)', 'FontSize', 20);
xlim([tday(100), tday(end)]);

% Exploratory plots of means and variances
figure;
subplot(4, 1, 1);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
plot(tday, ones(size(Rtrue)), '--', 'color', grey2, 'linewidth', 2);
ylabel('$R_{s}$', 'FontSize', 20);
grid off; box off; hold off;
xlim(xlim1);
subplot(4, 1, 2);
plot(tday, Fqm, '-', 'color', 'r', 'linewidth', 2);
hold on;
plot(tday, Fpm, '-', 'color', 'b', 'linewidth', 2);
ylabel('E[$1 - F_s(1)$]', 'FontSize', 20);
grid off; box off; hold off;
xlim(xlim1);
subplot(4, 1, 3);
plot(tday, Fqvar, '-', 'color', 'r', 'linewidth', 2);
hold on;
plot(tday, Fpvar, '-', 'color', 'b', 'linewidth', 2);
ylabel('var[$1 - F_s(1)$]', 'FontSize', 20);
grid off; box off; hold off;
xlim(xlim1);
subplot(4, 1, 4);
plot(tday, klm, '-', 'color', 'r', 'linewidth', 2);
ylabel('KL$(p | q)$', 'FontSize', 20);
grid off; box off;
xlabel('time, $s$ (days)', 'FontSize', 20);
xlim(xlim1);

%% Publication figures and data storage

% All change points as vertical lines
pltup = kron(tday(tup), ones(1, 2));
pltup = reshape(pltup, [2 nup]);
pltdwn = kron(tday(tdown), ones(1, 2));
pltdwn = reshape(pltdwn, [2 ndown]);

% Change-points to plot
figure('Renderer', 'painters', 'Position', [10 10 600 600]);
subplot(5, 1, [1 2]);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
plot(tday, ones(size(Rtrue)), '--', 'color', grey2, 'linewidth', 2);
plotCIRaw(tday', Rmf', Rlf', Ruf', 'b');
plotCIRaw(tday', Rm', Rl', Ru', 'r');
ylabel('$\hat{R}_s$', 'FontSize', 20);
grid off; box off; hold off;
xlim(xlim1);

subplot(5, 1, 3);
plot(tday, klm, '-', 'color', 'k', 'linewidth', 2);
hold on; xlim(xlim1); h = gca; ylim = h.YLim;
plot([pltup pltdwn], reshape(repmat(ylim, [1 nup+ndown]), [2 nup+ndown]), ...
    '--', 'color', grey2, 'linewidth', 2);
grid off; box off; hold off;
ylabel('D$(p_s | q_s)$', 'FontSize', 20);

subplot(5, 1, [4 5]);
plotCIRaw(tday', Fpm, Fpq(1, :)', Fpq(2, :)', 'b');
hold on;
plotCIRaw(tday', Fqm, Fqq(1, :)', Fqq(2, :)', 'r');
plot(tday, 0.5*ones(size(Rtrue)), '--', 'color', grey2, 'linewidth', 2);
xlim(xlim1); ylim = [0 1.1];
plot([pltup pltdwn], reshape(repmat(ylim, [1 nup+ndown]), [2 nup+ndown]), ...
    '--', 'color', grey2, 'linewidth', 2);
hold off; grid off; box off;
ylabel('P($R_s > 1$)', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['comb' namstr], 'fig');
    cd(thisDir);
end

% Change-points to plot
figure('Renderer', 'painters', 'Position', [10 10 600 600]);
subplot(5, 1, [1 2]);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on; xlim(xlim1); 
plot(tday, ones(size(Rtrue)), '--', 'color', grey2, 'linewidth', 2);
plotCIRaw(tday', Rmf', Rlf', Ruf', 'b');
plotCIRaw(tday', Rm', Rl', Ru', 'r'); 
h = gca; ylim = h.YLim;
plot([pltup pltdwn], reshape(repmat(ylim, [1 nup+ndown]), [2 nup+ndown]), ...
    '--', 'color', grey2, 'linewidth', 2);
ylabel('$\hat{R}_s$', 'FontSize', 20);
grid off; box off; hold off;
xlim(xlim1);

subplot(5, 1, 3);
plot(tday, klm, '-', 'color', 'k', 'linewidth', 2);
hold on; xlim(xlim1); h = gca; ylim = h.YLim;
plot([pltup pltdwn], reshape(repmat(ylim, [1 nup+ndown]), [2 nup+ndown]), ...
    '--', 'color', grey2, 'linewidth', 2);
grid off; box off; hold off;
ylabel('D$(p_s | q_s)$', 'FontSize', 20);

subplot(5, 1, [4 5]);
plot(tday, Fqm, '-', 'color', 'r', 'linewidth', 2);
hold on;
plot(tday, Fpm, '-', 'color', 'b', 'linewidth', 2);
plot(tday, 0.5*ones(size(Rtrue)), '--', 'color', grey2, 'linewidth', 2);
xlim(xlim1); ylim = [0 1.1];
plot([pltup pltdwn], reshape(repmat(ylim, [1 nup+ndown]), [2 nup+ndown]), ...
    '--', 'color', grey2, 'linewidth', 2);
hold off; grid off; box off;
ylabel('P($R_s > 1$)', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 20);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['comb2' namstr], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    clear('pRup', 'pstate', 'Rlowf', 'Rlow', 'Rhighf', 'Rhigh', 'js', 'kl',...
        'Rmean', 'Rmeanf', 'ttemp', 'tqtemp', 'tptemp', 'Lday', 'Iday');
    save(['proc' namstr '.mat']);
    cd(thisDir);
end