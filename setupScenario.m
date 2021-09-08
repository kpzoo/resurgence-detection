% Possible simulated epidemic scenarios
function simVals = setupScenario(scenNo)

% Specific scenario parameters
switch(scenNo)
    % Rs are distinct reprod nums, ts are switch points
    case 1
        % Rapidly controlled epidemic
        Rch = [2 0.5]; tch = 100;
    case 2
        % Square wave (recovered rapid control)
        Rch = [2.5 0.5 2.5]; tch = [70 230];
    case 3
        % Three stage control with sines
        %Rch = [4 0.6 2 0.2]; tch = [40 80 150];
        % For ebola
        Rch = [4 0.7 2 0.2]; tch = [40 150 220];
        % For covid
        %Rch = [4 0.7 2 0.4]; tch = [40 150 220];
    case 4
        % Exponential rise and fall
        tch = 30; Rch = [];
    case 5
        % Two stage control with noise
        Rch = [3 0.8 0.4]; tch = [15 180];
    case 6
        % Second wave dynamics (sines)
        Rch = [1.3 1.2]; tch = 3;
    case 7
        % Long period of low R between transmission and noise
        Rch = [2.5 0.1 1.3]; tch = [50 200];
    case 8
        % Exponential rise and fall and rise
        tch = [40 190]; Rch = [];
    case 9
        % Simple rise in transmission
        Rch = [1.2 2]; tch = 200;
end
simVals.Rch = Rch; simVals.tch = tch;