% Convert Rt to growth rates rt via Gamma
function rt = compGrowth(distvals, Rt)

% Assumptions and notes
% - gamma serial interval

if distvals.type == 2
    % Gamma distribution shape and scale
    shape = distvals.pm;
    scale = distvals.omega/shape;
    
    % Growth rate from Wallinga-Lipsitch
    rt = (Rt.^(1/shape) - 1)/scale;
end