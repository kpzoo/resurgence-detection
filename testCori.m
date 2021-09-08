% Test real time behaviour using EpiEstim maths
clearvars; clc; close all; tic;

% Assumptions and notes
% - complementary CDFs to determine is resurgence harder to detect                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     eraged and other state models handle noise

% Directory and where saving
thisDir = cd; saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);


%% Simple R = 1 plateau - how does P(R > 1) and P(R <= 1) behave

% Key variable lam = incidence 
plat = 1:100; 
p_plat = poisscdf(plat, plat);

figure;
plot(plat, p_plat, 'LineWidth', 2);
box off; grid off;
xlabel('epidemic size', 'FontSize', 20);
ylabel('P($R > 1$)', 'FontSize', 20);

%% Examine from Poisson statistics what is hardest to detect

% Test total infectiousness values and increments
lam = 1:50; del = 1:25; 
lenl = length(lam); lend = length(del);

% Resurgence Poiss CDF terms
Sup = zeros(lenl, lend);
for i = 1:lend
   Sup(:, i) = poisscdf(lam+del(i), lam);
end

figure;
plot(lam, Sup, 'LineWidth', 2);
box off; grid off;

% Test total infectiousness values and decrements
lam = 100:200; del = 20:100; 
lenl = length(lam); lend = length(del);

% Controlled Poiss CDF terms
Sdown = zeros(lenl, lend);
for i = 1:lend
   Sdown(:, i) = 1 - poisscdf(lam-del(i), lam);
end

figure;
plot(lam, Sdown, 'LineWidth', 2);
box off; grid off;

% Simple poisscdf test (from Cori)
l = 1:40; b = 1:10; lenb = length(b);
cori1 = zeros(length(l), lenb);
for i = 1:lenb
    cori1(:, i) = poisscdf(l+b(i),  l);
end


%% Simple inverse CDF visualisations

% Fix some quantiles of interest
pr = 0.1:0.01:0.9; lenr = length(pr);

% Epidemic sizes (mean)
lam = [10 20 50 100 150 200]; lenl = length(lam);
% Quantile values normalised by mean
Xr = zeros(lenr, lenl); 

% Get size perturbations necessary to achieve pr
for i = 1:lenr
    % This is analogous to R in mean for R > 1
    Xr(i, :) = poissinv(pr(i), lam)./lam;
end

figure;
stairs(pr', Xr(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
hold on;
stairs(pr', Xr(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(pr', Xr(:, end), 'LineWidth', 2, 'Color', 'r');
box off; grid off; hold off;
ylabel('relative size', 'FontSize', 20);
xlabel('P($R > 1$)', 'FontSize', 20);
xlim([pr(1) pr(end)]);

%% Inverse CDF but switch from R > 1 to R <= 1 for 0.5

% Size perturbations but with switch of meaning
pr = 0.01:0.01:0.99; lenr = length(pr);
% Epidemic sizes (mean)
lam = 20:20:200; lenl = length(lam);

% P(R > 1) and P(R <= 1) vectors
Ur = zeros(lenr, lenl); Dr = Ur;

% Get size perturbations necessary to achieve pr
for i = 1:lenr
    % This is analogous to R in mean for R > 1
    Ur(i, :) = poissinv(pr(i), lam)./lam - 1;
    Dr(i, :) = poissinv(1 - pr(i), lam)./lam - 1;
end

figure;
subplot(2, 1, 1);
stairs(pr', Ur(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
%stairs(pr', Ur(:, 2:end-1), '.', 'LineWidth', 2, 'Color', grey1, 'MarkerSize', 20);
hold on;
stairs(pr', Ur(:, 1), 'LineWidth', 2, 'Color', 'b');
%stairs(pr', Ur(:, 1), '.', 'LineWidth', 2, 'Color', 'b', 'MarkerSize', 20);
stairs(pr', Ur(:, end), 'LineWidth', 2, 'Color', 'r');
%stairs(pr', Ur(:, end), '.', 'LineWidth', 2, 'Color', 'r', 'MarkerSize', 20);
box off; grid off; hold off; set(gca, 'XGrid', 'on');
ylabel('relative size', 'FontSize', 20);
xlabel('P($R_s > 1 | I_1^s$)', 'FontSize', 20);
xlim([pr(1) pr(end)]);
subplot(2, 1, 2);
stairs(pr', Dr(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
hold on;
stairs(pr', Dr(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(pr', Dr(:, end), 'LineWidth', 2, 'Color', 'r');
box off; grid off; hold off; set(gca, 'XGrid', 'on');
ylabel('relative size', 'FontSize', 20);
xlabel('P($R_s \leq 1 | I_1^s$)', 'FontSize', 20);
xlim([pr(1) pr(end)]);

figure;
hold on;
for i = 2:lenl-1
    stairs(Ur(:, i), pr', 'LineWidth', 2, 'Color', grey1);
end
stairs(Ur(:, 1), pr', 'LineWidth', 2, 'Color', 'b');
stairs(Ur(:, end), pr', 'LineWidth', 2, 'Color', 'r');
h = gca; plot(h.XTick, 0.5*ones(size(h.XTick)), 'k--', 'LineWidth', 1);
plot(zeros(size(h.YTick)), h.YTick, 'k--', 'LineWidth', 1);
box off; grid off; hold off; set(gca, 'YGrid', 'on');
xlabel('relative epidemic size', 'FontSize', 20);
ylabel('P($R_s > 1 | I_1^s$)', 'FontSize', 20);
xlim([h.XTick(1) h.XTick(end)]);


%% Inverse CDF for given mean R = a plateau

% Key variable lam = incidence 
plat = 1:100; nplat = length(plat);
% Factors for Rmean
fac = 1.1:0.1:2; nfac = length(fac);

% For each Rmean contour get possible prob of plateau
p_plat = zeros(nfac, nplat);
for i = 1:nfac
    p_plat(i, :) = poisscdf(fac(i)*plat, plat);
end

figure;
stairs(plat', p_plat', 'LineWidth', 2);
box off; grid off;
xlabel('epidemic size', 'FontSize', 20);
ylabel('P($R > 1$)', 'FontSize', 20);


%% Plot of gamma distributions for various changes

% Changes as a % of total infectiousness
del = [0.1 0.25 0.5 0.75 1]; ndel = length(del);
% Total infectiousness + c term
lam = [25 100 400]; 
scale = 1./lam; nscale = length(lam);

% Shape (infecteds) for each lam and delta and mean
shape = zeros(ndel, nscale); meanGam = shape;
for i = 1:ndel
    shape(i, :) = (1 + del(i))*lam;
    meanGam(i, :) = gamstat(shape(i, :), scale);
end

% Set of R and gamma probabilities
Rset = 0.001:0.001:3; lenR = length(Rset);
Pgam = cell(ndel, nscale);
for i = 1:ndel
    for j = 1:nscale
        Pgam{i, j} = gampdf(Rset, shape(i, j), scale(j));
    end
end

% Colours for different delta
colGam = {'b', grey1, grey1, grey1, 'r'};

% Plot possible gamma posteriors
figure;
for i = 1:nscale
    subplot(nscale, 1, i);
    hold on;
    for j = 1:ndel
        plot(Rset, Pgam{j, i}, 'Color', colGam{j}, 'LineWidth', 2);
    end
    box off; grid off;
    ylabel('P($R_s | I_1^s$)', 'FontSize', 20);  
    title(['$\lambda_s$ = ' num2str(lam(i))]);
end
xlabel('$R_s$', 'FontSize', 20);   



%% Publishable figure

figure('Renderer', 'painters', 'Position', [10 10 600 500]);
for i = 1:nscale
    subplot(3, 2, 2*i-1);
    hold on;
    for j = 1:ndel
        plot(Rset, Pgam{j, i}, 'Color', colGam{j}, 'LineWidth', 2);
    end
    hold off; box off; grid off; set(gca, 'XGrid', 'on');
    ylabel('P($R_s | I_1^s$)', 'FontSize', 20);  
    title(['$\lambda_{\tau(s)}$ = ' num2str(lam(i))]);
end
xlabel('$R_s$', 'FontSize', 20); 

subplot(3, 2, [2 4 6]);
hold on;
for i = 2:lenl-1
    stairs(Ur(:, i), pr', 'LineWidth', 2, 'Color', grey1);
end
stairs(Ur(:, 1), pr', 'LineWidth', 2, 'Color', 'b');
stairs(Ur(:, end), pr', 'LineWidth', 2, 'Color', 'r');
h = gca; plot(h.XTick, 0.5*ones(size(h.XTick)), 'k--', 'LineWidth', 1);
plot(zeros(size(h.YTick)), h.YTick, 'k--', 'LineWidth', 1);
box off; grid off; hold off; set(gca, 'YGrid', 'on');
%xlabel('($i_{\tau(s)}\lambda_{\tau(s)}^{-1} - 1)$', 'FontSize', 20);
xlabel('$\Delta \lambda_{\tau(s)}$', 'FontSize', 20);
ylabel('P($R_s > 1 | I_1^s$)', 'FontSize', 20);
xlim([h.XTick(1) h.XTick(end)]);


