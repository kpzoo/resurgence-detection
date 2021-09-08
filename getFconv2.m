% Function to compute P(R > 1) under different weighted mean scenarios
function [F1, F2all, F3all, Fstruc] = getFconv2(w1, L1, del, Rgrid, dRval, f33)

% Assumptions and notes
% - Allows manipulation of incidence, total infectiousness and R in regions

% Delta to incidence on region 1, rest keep at R <= 1
I1 = (1 + del)*L1; ndel = length(del);

% Coefficients for Rs in other regions
f22 = 1; f23 = 1; 
% Total infectiousness weights
w23 = (1 - w1)/2; w33 = 1 - w1 - w23; 
% Total infectiousness
L22 = ((1-w1)/w1)*L1; L23 = (w23/w1)*L1; L33 = (w33/w1)*L1;

% Key P(R > 1) distributions
F1 = zeros(1, ndel); F2all = F1; F3all = F1;

% Distributions of region 1 and overall
for i = 1:ndel
    % Region 1 with and without weight
    p1w = gampdf(Rgrid, I1(i), w1/L1); p1w = p1w/sum(p1w); 
    p1 = gampdf(Rgrid, I1(i), 1/L1); p1 = p1/sum(p1);
    
    % Overall R distribution for p = 2, R2 = a22 (I2 = a2*L2)
    p22 = gampdf(Rgrid, f22*L22, (1-w1)/(L22)); p22 = p22/sum(p22);
    p2all = conv(p1w, p22)*dRval;
    p2all = interp1(2*Rgrid, p2all(1:2:end), Rgrid,'linear','extrap');
    p2all = p2all/(sum(p2all));
    
    % Overall R distribution for p = 3, R2 = a23, R3 = a33
    p23 = gampdf(Rgrid, f23*L23, w23/L23); p23 = p23/sum(p23);
    p3all = conv(p1w, p23)*dRval;
    p3all = interp1(2*Rgrid, p3all(1:2:end), Rgrid,'linear','extrap');
    p3all = p3all/(sum(p3all));

    % Normalise and repeat convolution
    p33 = gampdf(Rgrid, f33*L33, w33/L33); p33 = p33/sum(p33);
    p3all = conv(p3all, p33)*dRval;
    p3all = interp1(2*Rgrid, p3all(1:2:end), Rgrid,'linear','extrap'); 
    p3all = p3all/(sum(p3all));
    
    % P(R1 > 1) and overall P(R > 1)
    id = find(Rgrid == 1);
    F1(i) = 1 - sum(p1(1:id)); 
    F2all(i) = 1 - sum(p2all(1:id));
    F3all(i) = 1 - sum(p3all(1:id));
end

% Check consistency of weights
if ~(abs((1-w1) - L22/(L22 + L1)) < 10^-10 && abs(w23 - L23/(L23 + L1 + L33))...
        < 10^-10 && abs(w33 - L33/(L23 + L1 + L33)) < 10^-10)
    error('Weights are not sensible');
end

% Output structure
Fstruc.w2 = [w1, 1-w1];
Fstruc.w3 = [w1, w23, w33];
Fstruc.L2 = [L1 L22];
Fstruc.L3 = [L1 L23 L33];
Fstruc.f2 = f22;
Fstruc.f3 = [f23 f33];