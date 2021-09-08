
tday = 1:300; nday = length(tday);
L = 1.5; trise = 210; tfall = tday(end) - trise;
steep = logspace(log10(0.08), log10(0.4), 5); nst = length(steep);
Rrise = zeros(nst, nday); Rfall = Rrise;

for i = 1:nst
    Rrise(i, :) = 0.25 + L./(1 + exp(-steep(i)*(tday - trise)));
    Rfall(i, :) = 0.25 + L - L./(1 + exp(-steep(i)*(tday - tfall)));
end

figure;
plot(tday, Rrise + Rfall, 'LineWidth', 2);