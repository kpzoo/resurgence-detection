% Process real time delay data to get lag statistics
clearvars; clc; close all; tic;

% Assumptions and notes
% - logistic simulations to determine rise and fall delays in R
% - batch across 4 types of disease and 5 possible logistic functions

% Directory and where saving
thisDir = cd; saveFol = 'Data/logistic';
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Basic parameters for simulations

% Set number of epidemics
M = 1000; disp(['Simulating ' num2str(M) ' epidemics']);
% Define number of days to simulate
nday0 = 301; tday0 = 1:nday0;

% Possible infectious diseases
epiNo = 1:5; nEpi = length(epiNo);

% True R profiles
L = 1.5; trise = 210; tfall = tday0(end) - trise;
% Rising and falling coefficients
steep = logspace(log10(0.08), log10(0.4), 10);
nScen = length(steep); Rtrues = zeros(nScen, nday0);

% Logistic type functions that look like step recovery
for i = 1:nScen
    Rtrues(i, :) = 0.5 + L./(1 + exp(-steep(i)*(tday0 - trise))) +...
        L - L./(1 + exp(-steep(i)*(tday0 - tfall)));
end

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% Probability countours for checking F(1) and 1 - F(1)
pr = [0.5 0.75 0.95]; lenr = length(pr);
% Crossing times
tups = zeros(nEpi, nScen); tdowns = tups;

% Mean delays in crossing R > pr 
tdel_up = cell(nEpi, nScen); tdel_down = tdel_up;
% Serial interval distributions
Pserial = cell(1, nEpi);

%% Main simulation code for each disease with R estimation

for iii = 1:nEpi
    
    % For each disease consider a logistic scenario
    for ii = 1:nScen
        
        % Epidemics variables to store
        Iday = cell(1, M); Lday = Iday; tday = Iday;
        
        % Simulate M epidemic replicates
        for i = 1:M
            % Simulation parameters and warning
            Iwarn = 1;
            while Iwarn
                % Main epidemic simulation function
                [Iday{i}, Lday{i}, Rtrue, tday{i}, Iwarn, distvals, Pomega] =...
                    epiSimScenLogistic(nday0, epiNo(iii), Rtrues(ii, :));
                
                % Remove sequences with few cases at end
                idrem = sum(Iday{i}(end-5:end)) < 20;
                if idrem
                    Iwarn = 1;
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
        
        % Ensure single crossing
        nup = length(tup); ndown = length(tdown);
        if nup ~= 1 || ndown ~= 1
            error('Incorrect number of crossings');
        else
            % Store crossings
            tups(iii, ii) = tday(tup); tdowns(iii, ii) = tday(tdown);
            disp(['Downward crossing at: ' num2str(tday(tdown))]);
            disp(['Upward crossing at: ' num2str(tday(tup))]);
        end
        
        % R and incidence variables
        Rmean = cell(1, M); Rhigh = Rmean; Rlow = Rmean;
        Rmeanf = Rmean; Rhighf = Rmean; Rlowf = Rmean;
        % P(R > 1) for smoothing and filtering
        Fq = zeros(nday, M); Fp = Fq;
        
        for i = 1:M
            % Standard EpiFilter and EpiSmoother estimates
            [~, ~, ~, ~, pR, pRup, pstate] = ...
                runEpiFilter(Rgrid, m, eta, nday, p0, Lday{i}, Iday{i});
            [~, ~, ~, ~, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);
            
            % Complementary CDFs for R > 1
            id1 = find(Rgrid <= 1, 1, 'last');
            Fq(:, i) = 1 - sum(qR(:, 1:id1), 2)';
            Fp(:, i) = 1 - sum(pR(:, 1:id1), 2)';
        end
        
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
        
        % Store crossing information
        tdel_up{iii, ii} = tmeanup;  tdel_down{iii, ii} = tmeandown;
    end
    
    % Save serial interval
    Pserial{iii} = Pomega;
    % Progress updates
    disp(['Completed ' num2str(iii) ' of ' num2str(nEpi)]);
end

% Key data save
clear('tpset', 'tpset2', 'tqset', 'tqset2', 'prUp', 'pR', 'pstate', 'qR', 'Lday', 'Iday');
save(['batchLog_' num2str(nEpi) '_' num2str(nScen) '_' num2str(M)]);

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

%% Examine timing and plot

% Parse data at pr = 0.5 and 0.75
id1 = find(pr == 0.5); id2 = find(pr == 0.95);
up1 = zeros(nEpi, nScen); down1 = up1; up2 = up1; down2 = up1;

% Extract from overall cell matrix
for i = 1:nEpi
    for j = 1:nScen
        % Difference in time for two pr = P(R > 1) being achieved
        up1(i, j) = tdel_up{i, j}(id1); down1(i, j) = tdel_down{i, j}(id1);
        up2(i, j) = tdel_up{i, j}(id2); down2(i, j) = tdel_down{i, j}(id2);
    end   
end

% Colours for 5 diseases
cols = {'g', grey1, grey2, 'r', 'b'};

% Steepness against upward and downward differences at 0.5
figure;
subplot(2, 1, 1);
hold on;
for i = 1:nEpi
    plot(steep, up1(i, :), '.', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 30);
    plot(steep, down1(i, :), 'd', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', cols{i});
    %stairs(steep, up1(i, :), 'Color', cols{i}, 'LineWidth', 2);
    %stairs(steep, down1(i, :), 'Color', cols{i}, 'LineWidth', 2);
end
hold off; grid off; box off; set(gca, 'XGrid', 'on');
ylabel('$\Delta t_{50}$', 'FontSize', 20);
subplot(2, 1, 2);
hold on;
for i = 1:nEpi
    plot(steep, up2(i, :), '.', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 30);
    plot(steep, down2(i, :), 'd', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', cols{i});
end
hold off; grid off; box off; set(gca, 'XGrid', 'on');
ylabel('$\Delta t_{75}$', 'FontSize', 20);
xlabel('steepness, $k$', 'FontSize', 20);

%% Publishable

% Four panel figure of all results
figure('Renderer', 'painters', 'Position', [10 10 600 600]);
% True R numbers
subplot(3, 2, 1);
plot(tday0, Rtrues(2:end, :), 'Color', grey1, 'LineWidth', 2);
hold on;
plot(tday0, Rtrues(1, :), 'Color', 'r', 'LineWidth', 2);
plot(tday0, Rtrues(end, :), 'Color', 'b', 'LineWidth', 2);
plot(tday0, ones(size(tday0)), 'k--', 'LineWidth', 1);
hold off; grid off; box off;
ylabel('$R_s$', 'FontSize', 20); ylim([0 2.5]);
xlabel('time, $s$ (days)', 'FontSize', 20);

% Serial intervals
subplot(3, 2, 2);
hold on;
for i = 1:nEpi
    plot(tday0, Pserial{i}, 'Color', cols{i}, 'LineWidth', 2);
end
hold off; grid off; box off;
ylabel('$w_u$', 'FontSize', 20); xlim([0 45]);
xlabel('time, $u$ (days)', 'FontSize', 20);

% P(R > 1) at 50%
subplot(3, 2, [3 5]);
hold on;
for i = 1:nEpi
    plot(steep, up1(i, :), '.', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 25);
    plot(steep, down1(i, :), 'd', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', cols{i});
end
hold off; grid off; box off; xlim([0.07 0.41]);
ylabel('$\Delta t_{50}$', 'FontSize', 20);
xlabel('steepness, $k$', 'FontSize', 20);

% P(R > 1) at 95%
subplot(3, 2, [4 6]);
hold on;
for i = 1:nEpi
    plot(steep, up2(i, :), '.', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 25);
    plot(steep, down2(i, :), 'd', 'Color', cols{i}, 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', cols{i});
end
hold off; grid off; box off; xlim([0.07 0.41]);
ylabel('$\Delta t_{95}$', 'FontSize', 20);
xlabel('steepness, $k$', 'FontSize', 20);



