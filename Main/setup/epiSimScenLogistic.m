% Simulate epidemics via renewal model for different diseases
function [Iday, Lam, Rtrue, tday, Iwarn, distvals, Pomega] = epiSimScenLogistic(nday, epiNo, Rtrue)

% Assumptions and notes
% - includes different disease type scenarios, called by logisticSi.m
% - removes 'burn-in' of first 20 days, epidemic size < 100
% - various R trajectories pre-specified via logistic functions

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
