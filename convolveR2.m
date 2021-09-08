% Plot convolution of Rj distributions
clearvars; clc; close all; tic;

% Assumptions and notes
% - assumes each Rj has same component distribution
% - for now only including gamma distributions
% - region 1 always has the resurgence

% Directory of some main code and plotting options
thisDir = cd; cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Two gamma distributions 

% Support of all distributions
dR = 0.001; Rgrid = 0.001:dR:6; lenR = length(Rgrid);

% Total infectiousness in each region scaled by a
a = 2; L1 = 20; L2 = a*L1; I2 = L2;
% Weights for the regions based on a
w1 = L1/(L1 + L2); w2 = 1 - w1;

% Delta to incidence on region 1, rest keep at R = 1
del = [0.1 0.25 0.5 0.75 1]; 
I1 = (1 + del)*L1; ndel = length(del);

% Represent each region R distribution with a gamma pdf
p1 = zeros(ndel, lenR); pall = p1; 
% Distribution of region 2 and values
p2dist = @(x)gampdf(x, I2, w2/L2); p2 = p2dist(Rgrid)/sum(p2dist(Rgrid));

% Statistics for testing
m1 = zeros(1, ndel); v1 = m1; mall = m1; vall = m1;
% P(R > 1) quantities
F1 = zeros(1, ndel); Fall = F1;

% Distributions of region 1 and overall
for i = 1:ndel
    % Region 1 with and without weight
    p1w = @(x)gampdf(x, I1(i), w1/L1); p1dist = @(x)gampdf(x, I1(i), 1/L1);
    p1(i, :) = p1dist(Rgrid); p1(i, :) = p1(i, :)/sum(p1(i, :));
    
    % Overall R distribution is convolution
    pdist = conv(p1w(Rgrid)/sum(p1w(Rgrid)), p2)*dR;
    % Convolution is twice as long so shorten
    pdist = interp1(2*Rgrid, pdist(1:2:end), Rgrid,'linear','extrap'); 
    pall(i, :) = pdist/sum(pdist);
    
    % Statistics of region 1 and overall
    [m1(i), v1(i)] = gamstat(I1(i), 1/L1);
    mall(i) = pall(i, :)*Rgrid'; vall(i) = pall(i, :)*(Rgrid.^2)' - mall(i)^2;
    
    % P(R1 > 1) and overall P(R > 1)
    id = find(Rgrid == 1);
    F1(i) = 1 - sum(p1(i, 1:id)); Fall(i) = 1 - sum(pall(i, 1:id));
end

% Save key variables
p1_2 = p1; pall_2 = pall; F1_2 = F1; Fall_2 = Fall; w1_2reg = w1;
m1_2 = m1; v1_2 = v1; mall_2 = mall; vall_2 = vall; a_2 = a;

%% Three gamma distributions (but same w1)

% Total infectiousness in each region scaled by a
a = [1 1]; L1 = 20; L2 = a(1)*L1; L3 = a(2)*L1; I2 = L2; I3 = 0.8*L3;
% Weights for the regions based on a
w1 = L1/(L1 + L2 + L3); w2 = L2/(L1 + L2 + L3); w3 = 1 - (w1 + w2);

% Delta to incidence on region 1, rest keep at R <= 1
del = [0.1 0.25 0.5 0.75 1]; 
I1 = (1 + del)*L1; ndel = length(del);

% Represent each region R distribution with a gamma pdf
p1 = zeros(ndel, lenR); pall = p1; 
% Distribution of region 2, 3 and values
p2dist = @(x)gampdf(x, I2, w2/L2); p2 = p2dist(Rgrid)/sum(p2dist(Rgrid));
p3dist = @(x)gampdf(x, I3, w3/L3); p3 = p3dist(Rgrid)/sum(p3dist(Rgrid));

% Statistics for testing
m1 = zeros(1, ndel); v1 = m1; mall = m1; vall = m1;
% P(R > 1) quantities
F1 = zeros(1, ndel); Fall = F1;

% Distributions of region 1 and overall
for i = 1:ndel
    % Region 1 with and without weight
    p1w = @(x)gampdf(x, I1(i), w1/L1); p1dist = @(x)gampdf(x, I1(i), 1/L1);
    p1(i, :) = p1dist(Rgrid); p1(i, :) = p1(i, :)/sum(p1(i, :));
    
    % Overall R distribution is convolution of all regions
    pdist = conv(p1w(Rgrid)/sum(p1w(Rgrid)), p2)*dR;
    % Convolution is twice as long so shorten
    pdist = interp1(2*Rgrid, pdist(1:2:end), Rgrid,'linear','extrap'); 
    % Normalise and repeat convolution
    pdist = pdist/sum(pdist); pdist = conv(pdist, p3)*dR;
    pdist = interp1(2*Rgrid, pdist(1:2:end), Rgrid,'linear','extrap'); 
    
    % Total R distribution for overall regions
    pall(i, :) = pdist/sum(pdist);
    
    % Statistics of region 1 and overall
    [m1(i), v1(i)] = gamstat(I1(i), 1/L1);
    mall(i) = pall(i, :)*Rgrid'; vall(i) = pall(i, :)*(Rgrid.^2)' - mall(i)^2;
    
    % P(R1 > 1) and overall P(R > 1)
    id = find(Rgrid == 1);
    F1(i) = 1 - sum(p1(i, 1:id)); Fall(i) = 1 - sum(pall(i, 1:id));
end

% Save key variables
p1_3 = p1; pall_3 = pall; F1_3 = F1; Fall_3 = Fall; w1_3reg = w1;
m1_3 = m1; v1_3 = v1; mall_3 = mall; vall_3 = vall; a_3 = a;


%% Simulate P(R1 > 1) and P(R > 1) for three key scenarios

% Delta to incidence on region 1, rest keep at R <= 1
del = -0.95:0.05:0.95; ndel = length(del);

% Case 1: multiple w1 values, fixed Lj and Rj
w1_a = logspace(log10(1/20), log10(19/20), 9); len1 = length(w1_a); L1_a = 20; 
% Store results of CDFs
F1_a = zeros(len1, ndel); F2all_a = F1_a;
% Obtain CDF values
for i = 1:len1
    [F1_a(i, :), F2all_a(i, :), ~, ~] = getFconv2(w1_a(i), L1_a, del, Rgrid, dR, 1);
end

% Case 2: change L1 values, fixed Rj and w1
w1_b = 1/3; L1_b = logspace(log10(20), log10(80), 20); len2 = length(L1_b);  
% Store results of CDFs
F1_b = zeros(len2, ndel); F2all_b = F1_b;
% Obtain CDF values
for i = 1:len2
    [F1_b(i, :), F2all_b(i, :), ~, ~] = getFconv2(w1_b, L1_b(i), del, Rgrid, dR, 1);
end

% Case 3: fix w1, L1 but varying R3 for 3rd region
L1_c = 20; w1_c = 1/3; f3_c = linspace(0.5, 1.2, 8); len3 = length(f3_c);

% Store results of CDFs
F1_c = zeros(len3, ndel); F3all_c = F1_c;
% Obtain CDF values
for i = 1:len3
    [F1_c(i, :), ~, F3all_c(i, :), ~] = getFconv2(w1_c, L1_c, del, Rgrid, dR, f3_c(i));
end

% For plotting try quantiles for middle right panel
Q1b = quantile(F1_b, [0.025, 0.5, 0.975]);
Qallb = quantile(F2all_b, [0.025, 0.5, 0.975]);

%% Publishable figure

% First panel on convolutions with same p1
figure('Renderer', 'painters', 'Position', [10 10 600 500]);
subplot(3, 2, 1);
plot(Rgrid, p1_3(2:end-1, :), 'Color', grey1, 'LineWidth', 2);
hold on;
plot(Rgrid, p1_3(1, :), 'Color', 'b', 'LineWidth', 2);
plot(Rgrid, p1_3(end, :), 'Color', 'r', 'LineWidth', 2);
hold off; box off; grid off; set(gca, 'XGrid', 'on');
xlim([0 3]); ylabel('P($R_s[1] | I_1^s$)', 'FontSize', 20);  
subplot(3, 2, 3);
plot(Rgrid, pall_2(2:end-1, :), 'Color', grey1, 'LineWidth', 2);
hold on; 
xlim([0 3]); ylabel('P($\bar{R}_s | p = 2$)', 'FontSize', 20);  
plot(Rgrid, pall_2(1, :), 'Color', 'b', 'LineWidth', 2);
plot(Rgrid, pall_2(end, :), 'Color', 'r', 'LineWidth', 2);
hold off; box off; grid off; set(gca, 'XGrid', 'on');
subplot(3, 2, 5);
plot(Rgrid, pall_3(2:end-1, :), 'Color', grey1, 'LineWidth', 2);
hold on;
plot(Rgrid, pall_3(1, :), 'Color', 'b', 'LineWidth', 2);
plot(Rgrid, pall_3(end, :), 'Color', 'r', 'LineWidth', 2);
hold off; box off; grid off; set(gca, 'XGrid', 'on');
xlim([0 3]); ylabel('P($\bar{R}_s | p = 3$)', 'FontSize', 20);  
xlabel('$R_s$', 'FontSize', 20); 

% Second panel are P(R > 1) curves
subplot(3, 2, 2);
plot(del, F2all_a(2:end, :), 'Color', grey1, 'LineWidth', 2);
hold on;
plot(del, F2all_a(1, :), 'Color', 'b', 'LineWidth', 2);
plot(del, F1_a(1, :), 'Color', 'r', 'LineWidth', 2);
h = gca; plot(h.XTick, 0.5*ones(size(h.XTick)), 'k--', 'LineWidth', 1);
plot(zeros(size(h.YTick)), h.YTick, 'k--', 'LineWidth', 1);
hold off; box off; grid off; %set(gca, 'YGrid', 'on');
ylabel('P($\bar{R}_s > 1 | \alpha_1$)', 'FontSize', 20);  
ylim([0 1.02]);
subplot(3, 2, 4);
% plot(del, F1_b(2:end-1, :), '.', 'Color', grey1, 'LineWidth', 2, 'MarkerSize', 10);
% hold on;
% plot(del, F2all_b(2:end-1, :), '--', 'Color', grey1, 'LineWidth', 2);
% plot(del, F1_b(1, :), '.', 'Color', 'b', 'LineWidth', 2, 'MarkerSize', 10);
% plot(del, F1_b(end, :), '.', 'Color', 'r', 'LineWidth', 2, 'MarkerSize', 10);
% plot(del, F2all_b(1, :), '--', 'Color', 'b', 'LineWidth', 2);
% plot(del, F2all_b(end, :), '--', 'Color', 'r', 'LineWidth', 2);

plotCIRaw(del', Q1b(2, :)', Q1b(1, :)', Q1b(3, :)', 'b');
hold on;
plotCIRaw(del', Qallb(2, :)', Qallb(1, :)', Qallb(3, :)', 'r');

h = gca; plot(h.XTick, 0.5*ones(size(h.XTick)), 'k--', 'LineWidth', 1);
plot(zeros(size(h.YTick)), h.YTick, 'k--', 'LineWidth', 1);
hold off; box off; grid off; %set(gca, 'YGrid', 'on');
ylabel('P($\bar{R}_s > 1 | \lambda_{\tau(s)}$)', 'FontSize', 20);
ylim([0 1.02]);
subplot(3, 2, 6);
plot(del, F3all_c(2:end-1, :), 'Color', grey1, 'LineWidth', 2);
hold on;
plot(del, F3all_c(1, :), 'Color', 'b', 'LineWidth', 2);
plot(del, F3all_c(end, :), 'Color', 'r', 'LineWidth', 2);
plot(del, F1_c(1, :), 'k', 'LineWidth', 2);
h = gca; plot(h.XTick, 0.5*ones(size(h.XTick)), 'k--', 'LineWidth', 1);
plot(zeros(size(h.YTick)), h.YTick, 'k--', 'LineWidth', 1);
hold off; box off; grid off; %set(gca, 'YGrid', 'on');
ylabel('P($\bar{R}_s > 1 | R_s[3]$)', 'FontSize', 20);
xlabel('$\Delta \lambda_{\tau(s)}[1]$', 'FontSize', 20); 
ylim([0 1.02]);