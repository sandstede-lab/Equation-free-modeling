%% Generates data for traffic
% genTrafficData(5000, 0.96, 1.1, 0, 4.5, 200, 500, '../data/data_full.csv',  '../data/data_headways.csv')
% dataPoints-   number of different simulations to run
% v0Min-        min v0 to run each simulation
% v0Max-        max v0 to run each simulation
% muMin-        min perturbation for sin initial condition
% muMax-        max perturnbation for sin initial condition
%               actual mu will be uniformly distributed between muMin and muMax
% tMin =        the minimum time to record the simulation at
% tMean-        mean of the exponential distribution for the times
% dataFile-     output file name for full data, should be .csv
% hwaysFile-    output file name for headways, should be .csv
function genTrafficData(dataPoints, v0Min, v0Max, muMin, muMax,tMin, tMean, dataFile, hwaysFile)

options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
rng(18); % set random seed
h = 2.4;
len = 60;
numCars = 30;
cars = zeros(numCars*2, 1);
allData = zeros(dataPoints, 2*numCars + 3);

for j=1:dataPoints
    if mod(j, 100) ==0
        disp(j);
    end
    
    % randomly generate initial condition, v0, t
    t= exprnd(tMean)+tMin; % drawn from shifted exponential
    v0 = (v0Max-v0Min)*rand()+v0Min;
    mu = muMin + rand()*(muMax - muMin);    

    % initialize the cars
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(h, len/numCars, v0);
    end

    % evolve cars
    [~,allTime] = ode45(@microsystem,[0 t],cars, options, [v0 len h]);  
    
    % record the evolution and parameters
    allData(j, 1:2*numCars) = allTime(end,:);       
    allData(j, end) = v0;
    allData(j, end-1) = mu;
    allData(j, end-2) = t;

    % save the snapshots of the simulation
    if(mod(j,200)==0)
        writematrix(allData, dataFile);  
    end
end

% save the headways separately 
hways = getHeadways(allData(:,1:numCars)', len);
writematrix(allData, dataFile);
writematrix(hways, hwaysFile);

end
