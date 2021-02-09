%% Generates data for traffic
% dataPoints-   number of different simulations to run
% v0Min-        min v0 to run each simulation
% v0Max-        max v0 to run each simulation
% muMin-        min perturbation for sin initial condition
% muMax-        max perturnbation for sin initial condition
%               actual mu will be uniformly distributed between muMin and muMax
% tMin =        the minimum time to record the simulation at
% tMean-        mean of the exponential distribution for the times
% dataFile-     output file name, should be .mat
function genTrafficData(dataPoints, v0Min, v0Max, muMin, muMax,tMin, tMean, dataFile)

options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
h = 2.4;
len = 60;
numCars = 30;
cars = zeros(numCars*2, 1);
data = zeros(dataPoints, 2*numCars + 3);

for j=1:dataPoints
    disp(j);
    % randomly generate initial condition, v0, t
    t= exprnd(tMean)+tMin; % drawn from shifted exponential
    v0 = (v0Max-v0Min)*rand()+v0Min;
    mu = muMin + rand()*(muMax - muMin);    

    % initialize the cars
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(h, len/numCars, v0);
    end

    [~,allTime] = ode45(@microsystem,[0 t],cars, options, [v0 len h]);           % evolve cars
    data(j, 1:2*numCars) = allTime(end,:);                                       % record the evolution
    data(j, end) = v0;
    data(j, end-1) = mu;
    data(j, end-2) = t;

    if(mod(j,20)==0)
        save(dataFile,'data');  % save the snapshots of the simulation
    end
end


save(dataFile,'data');  % save the snapshots of the simulation

end
