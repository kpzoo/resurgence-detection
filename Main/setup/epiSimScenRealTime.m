% Simulate epidemics via renewal model for different diseases
function [Iday, Lam, Rtrue, tday, Iwarn, distvals] = epiSimScenRealTime(scenNo, nday, epiNo, simVals)

% Assumptions and notes
% - includes different disease type scenarios
% - removes 'burn-in' of first 20 days, epidemic size < 100
% - various R trajectories but only gamma serial intervals

% Possible scenarios available - must match calling function
scenNam = {'control', 'square-wave', 'cascade', 'boom-bust', 'filtered', 'waves', 'noise valley', 'boom-bust-boom', 'rising'};
disp(['True R scenario: ' scenNam{scenNo}]);
epiNam = {'Marburg', 'MERS', 'Measles', 'COVID-19', 'EVD'};
disp(['Disease: ' epiNam{epiNo}]);

%% Serial interval specification

% Hyerparameters of serial distribution (must be gamma)
switch(epiNo)
    case 1
        % Marburg distribution from van Kerkhove 2015
        omega = 9; pm = (omega^2)/(5.4^2); 
    case 2
        % MERS distribution from Cauchemez 2016
        omega = 6.8; pm = (omega^2)/(4.1^2);
    case 3
        % Measles distribution from Cori 2013
        omega = 14.9; pm = (omega^2)/(3.9^2);
    case 4
        % COVID-19 distribution from Ferguson 2020
        omega = 6.5; pm = (1/0.65)^2;
    case 5
        % EVD distribution from van Kerkhove 2015
        omega = 15.3; pm = (omega^2)/(9.3^2);
end
distvals.omega = omega; distvals.pm = pm; distvals.type = 2;

% Serial distribution over all tday
serial = serialDistrTypes(nday, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);


%% True R number scenarios to simulate

% Parameters for changepoints in R trajectory
Rch = simVals.Rch; tch = simVals.tch;
% Variable for true R
Rtrue = zeros(1, nday);

% Functions for scenarios: R on a daily basis
switch(scenNo)
     case 1
        % Rapidly controlled epidemic
        Rtrue(1:tch) = Rch(1);
        Rtrue(tch+1:end) = Rch(2);
    case 2
        % Rapid control that recovers
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
    case 3
        % Three stage control with fluctuations
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:tch(3)) = Rch(3);
        Rtrue(tch(3)+1:end) = Rch(4);
        Rtrue = Rtrue + 0.3*cosd(2*(1:nday));
        Rtrue(Rtrue <= 0) = 0.1;
    case 4
        % Exponential rise and fall
        trise = 1:tch; tfall = tch+1:nday;
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.02*(1:tch)); Rmax = Rtrue(tch);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.008*(tfall - tch));
    case 5
        % Two stage control with filtered noise
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
        % Add noise and filter
        Rtrue = Rtrue + normrnd(0.5, 2, [1 nday]);
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
    case 6
        % Second (sine) wave dynamics
        Rtrue = Rch(1) + Rch(2)*sind(tch(1)*(1:nday));
    case 7
        % Long period of low R between transmission and noise
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
        % Add noise and filter
        Rtrue = Rtrue + normrnd(0.5, 2, [1 nday]);
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
        % Check for noise
        if any(Rtrue < 0)
            Rtrue(Rtrue < 0) = 0;
        end
    case 8
        % Exponential rise and fall then rise
        tch1 = tch(1); tch2 = tch(2);
        trise = 1:tch1; tfall = tch1+1:tch2; 
        triseSec = tch2+1:nday; % second wave
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.03*(1:tch1)); Rmax = Rtrue(tch1);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.015*(tfall - tch1));
        % Second wave
        Rtrue(triseSec) = Rtrue(tch2)*exp(0.02*(triseSec - tch2));
    case 9
        % Simple increase in transmissibility all at R > 1
        Rtrue(1:tch) = Rch(1);
        Rtrue(tch+1:end) = Rch(2);
end


%% Simulate epidemic based on chosen parameters

% Daily incidence and infectiousness
Iday = zeros(1, nday); Lam = Iday; 
% Initialise epidemic and warning
Iday(1) = 10; Iwarn = 0; 

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness
    Lam(i) = sum(Iday(i-1:-1:1).*Pomegat);    
    % Renewal incidence
    Iday(i) = poissrnd(Rtrue(i)*Lam(i));
end

% Remove start-up 20 days
idz = 20:nday; tday = idz;
% Adjusted vectors - including tday
Iday = Iday(idz); Rtrue = Rtrue(idz); Lam = Lam(idz);

% Remove small epidemics
if sum(Iday) < 500
    Iwarn = 1;
    disp(['Sum is ' num2str(sum(Iday))]);
end

% Add simple offset of zeros
idoff = 1:(length(Iday) - simVals.offset);
Iday = [zeros(1, simVals.offset) Iday(idoff)];
Rtrue = [zeros(1, simVals.offset) Rtrue(idoff)];
Lam = [zeros(1, simVals.offset) Lam(idoff)];
