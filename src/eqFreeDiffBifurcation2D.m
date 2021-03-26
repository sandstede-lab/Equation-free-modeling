function eqFreeDiffBifurcation2D()
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 30;           % number of cars
full = true;
tskip = 300;
nPeriod = 10;
k = 8;
numEigvecs = 2;
weight = 5;
foptions = optimoptions(@fsolve, 'Display', 'iter');

%% load diffusion map data
if full
    hways = readmatrix('../data/data_headways.csv');
    diffMap2D = DiffusionMap(hways, numEigvecs, weight);
    stepSize = 0.001; 
    steps = 96;
else
    newData = readmatrix('../data/1000data2D.csv');
    diffMap2D = DiffusionMap(newData, numEigvecs, weight);
    stepSize = 0.002; 
    steps = 182;
end

% load the reference states and previous bifurcation diagram
bif = readmatrix('../results/microBif.csv');
embed = diffMap2D.restrict(bif(1:numCars,:));
vel = bif(end, :);
n = sqrt(embed(1,:).^2 + embed(2,:).^2);
T = -numCars./bif(numCars+ 1,:);
% save embedding of micro system
if full
    writematrix([n', T', vel'], '../results/micro2D.csv');
else
    writematrix([n', T', vel'], '../results/1000micro2D.csv');
end

start = length(vel);
change = -1;
v0_base2 = vel(start + change);
v0_base1 = vel(start);
T1 = T(start);
T2 = T(start + change);
p2 = norm(embed(:, start + change));
p1 = norm(embed(:, start));
coord2 = [p2, T2, v0_base2];
coord1 = [p1, T1, v0_base1];
rayAngle = 0;

% draw the bifurcation diagrams
figure;

subplot(1,2,1);hold on;
scatter(vel, n, 50, 'r.');
xlabel('v_0');
ylabel('\rho');
scatter(v0_base1, p1, 400,'k.'); drawnow;
scatter(v0_base2, p2, 400,'k.'); drawnow;

subplot(1,2,2); hold on;
scatter(n, T, 50, 'r.');
scatter(p1,T1, 400, 'k.'); drawnow;
scatter(p2, T2, 400, 'k.'); drawnow;
xlabel('\rho');
ylabel('T');

%% initialize secant continuation
bif = zeros(4,steps);

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    
    w = coord2 - coord1;
    scaledW = w ./ coord2;
    newGuess = coord2 + stepSize *(w/norm(w)); % first guess on the secant line
    scaledNewGuess = newGuess ./ coord2;
    
    %% Newton's method
    u = fsolve(@(u)FW(u, diffMap2D, rayAngle, coord2, scaledW, scaledNewGuess), newGuess, foptions);
    [p, sig] = ler( [cos(rayAngle)*u(1); sin(rayAngle)*u(1)], diffMap2D, u(3));
    bif(1, iEq) = norm(p);
    bif(2:3,iEq) = u(2:3)';                                            % save the new solution
    
    %% reset the values for the arc length continuation
    coord1 = coord2;
    coord2 = u;
    bif(4,iEq) = sig; % save the standard deviation of the headways to compare
    
    if full
        writematrix(bif, '../results/bifurcation2D.csv');
    else
        writematrix(bif, '../results/1000bifurcation2D.csv');
    end
    
    % plot the diagram
    subplot(1,2,1); hold on;
    scatter(bif(3, iEq), bif(1,iEq), 400, 'b.'); drawnow;
    subplot(1,2,2); hold on;
    scatter(bif(1,iEq), bif(2,iEq), 400, 'b.'); drawnow;
    
end

%% function to zero
% PARAMETERS:
% u         - the pair of radius and period (sigma,T,v0) that we're trying to find with
% diffMap2D - diffusion map object
% coord2    - secant coordinate for normalization
% ray       - angle on the embedding
% scaledW   - scaled secant direction
% scaledNewGuess - scaled original prediction
% RETURNS:
% out  - the distance between the input point and the point determined by
%               lifting, evolving for T time, and restricting and dot
%               product with the secant line
% onplane   - whether this is on the plane orthogonal to w at the starting guess
    function out = periodicDistance(u, v, diffMap2D, ray)
        
        inputPoint = [cos(ray)*u(1); sin(ray)*u(1)];
        timeEvolve = nPeriod*u(2);
        [r0, ~, r1] = ler(inputPoint, diffMap2D, v, timeEvolve);
        
        % convert to polar coordinates
        angle0 = atan2(r0(2), r0(1));
        angle1 = atan2(r1(2), r1(1));
        out = [norm(r1) - norm(r0), angle1-angle0];
    end

    function F = FW(u, diffMap2D, ray, coord2, scaledW, scaledNewGuess)
        F = zeros(1,3);
        F(1:2) = periodicDistance(u(1:2), u(3), diffMap2D, ray);
        scaledU = u ./ coord2;
        F(3) = dot(scaledW, scaledU - scaledNewGuess);
    end
%% lift, evolve, restrict
% PARAMETERS:
% newval     - the current value of the macrovariables,
%                 used to seed the lifting
% diffMap2D  - the diffuson map object
% t          - the duration to evaluate the lifted parameters
% v0         - the optimal velocity parameter for this state
% RETURNS:
% sigma      - the macro-state of the headways after evolving for t
    function [coord, sigma, coord2] = ler(newval, diffMap2D, v0, t)
        options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
        
        % healing step
        liftedHeadway = diffMap2D.lift(newval, k);
        lifted = [hwayToPos(liftedHeadway) ; optimalVelocity(h, liftedHeadway, v0)];
        [~,evo] = ode45(@microsystem,[0 tskip],lifted, options,[v0 len h]);
        evoCars = getHeadways(evo(end,1:numCars)',len);
        coord = diffMap2D.restrict(evoCars);
        sigma = std(evoCars(1:numCars));
        
        if nargin > 3
            [~,evo2] = ode45(@microsystem,[0 t],evo(end,:), options,[v0 len h]);
            evoCars2 = getHeadways(evo2(end, 1:numCars)',len);
            sigma = std(evoCars2(1:numCars));
            coord2 = diffMap2D.restrict(evoCars2);
        end
    end

end
