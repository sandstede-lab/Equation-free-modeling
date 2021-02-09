h = 2.4;
len = 60;
numCars = 30;

v0_base1 = 1.14;
v0_base2 = 1.135;

%{
options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
%% initialize car positions and velocities
cars_1 = zeros(2*numCars, 1);
cars_2 = zeros(2*numCars, 1);
mu = .1;
finalTime  = 50000;
for i = 1:numCars
    cars_1(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
    cars_1(i+numCars) = optimalVelocity(len/numCars, v0_base1);
    
    cars_2(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
    cars_2(i+numCars) = optimalVelocity(len/numCars, v0_base2);
end
%% save new reference state (if you want)
tic;
[~,allTime_1] = ode45(@(t,y)microsystem(t,y, [v0_base1 len h]),[0 finalTime],cars_1, options);
ref_1 = allTime_1(end,:)';
[~,allTime_2] = ode45(@(t,y)microsystem(t,y, [v0_base2 len h]),[0 finalTime],cars_2, options);
ref_2 = allTime_2(end,:)';
toc;
save('refStats.mat','ref_1','ref_2');
%}

%% or you can load them instead
load('refStats.mat','ref_1','ref_2');

hways1 = getHeadways(ref_1(1:numCars));
hways2 = getHeadways(ref_2(1:numCars));

%% Bifurcation on micro
steps = 200;
stepSize = 0.025;
cguess = -0.8581;

foptions = optimoptions(@fsolve, 'TolFun',1e-15,'TolX',1e-15, 'Display', 'off', ...
    'SpecifyObjectiveGradient', true, 'CheckGradient', true, ...
    'FiniteDifferenceType', 'central');

[~, max1] = max(hways1);
[~, max2] = max(hways2);
hways1 = circshift(hways1, [-max1  + 10, 0]);    % shift so max headway is in the middle-ish
hways2 =circshift (hways2, [-max2 + 10, 0]);
discreteScaling = 1;

% add more discretization points for smoother curve
hways1 = interp1(1:numCars,hways1,...
    linspace(1,numCars,numCars * discreteScaling))';
hways2 = interp1(1:numCars,hways2,...
    linspace(1,numCars,numCars * discreteScaling))';

sys1 = [hways1; cguess; 0; v0_base1];
sys2 = [hways2; cguess; 0; v0_base2];

bif  = zeros(numCars + 3, steps);

figure; hold on;
scatter(v0_base1, std(hways1), 100,'r*');
scatter(v0_base2, std(hways2), 100, 'r*');

for iMic = 1:steps
    w = sys2 - sys1;        % tangent vector
    change =  stepSize*(w/norm(w));
    newGuess = sys2 + change;       % initial guess for newton solver
    lastGuess = sys2;               % reference state for phase condition
    [u , ~] = fsolve(@(sys2)FW(sys2, lastGuess, numCars, len, 1/1.7, w, newGuess),...
        newGuess, foptions);
    r = u(1:discreteScaling: end  - 3);
    bif(:,iMic) = [r; u(end - 2: end)];
    
    % set variables for next iteration
    sys1 = sys2;
    sys2 = u;
    
    scatter(bif(end,iMic), std(r), 400, 'b.'); drawnow;
end

%save('microBif.mat', 'bif');
%% plot results
figure;
scatter(bif(end,:), std(bif(1:numCars,:)),50,'b.');

%% function to minimize with fsolve
% var       - state to vary in order to minimize fw
% ref       - reference state for phase condition
% nCars     - number of cars
% len       - track length
% tau       - momentum constant
% W         - tangent vector
% initGuess - first guess (point on tangent vector)
    function [fw,J] = FW(var, ref, nCars, len, tau, W, initGuess)
        v0 = var(end);
        fw = zeros(nCars*discreteScaling+3,1);
        [fw(1:nCars*discreteScaling+2), J1] = microFJ(var(1:end-1), ref, nCars,...
            discreteScaling*nCars , len, v0, tau);
        fw(end) = W' * (var-initGuess);
        J2 = W';
        J = [J1 ; J2]; 
    end
