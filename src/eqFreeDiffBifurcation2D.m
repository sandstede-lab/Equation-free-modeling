% eqFreeDiffBifurcation2D(true)
% computes the bifrucation diagram for the traffic system using original
% headways and a two-dimensional diffusion map
%
% full  - if true uses the 5000 point diffusion map
%       - if false uses the 1000 point diffusion map
function eqFreeDiffBifurcation2D(full)
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 30;           % number of cars
tskip = 300;            % tiime to evolve to slow manifold
nPeriod = 7;            % time to evolve for finite difference
k = 8;                  % number of points to use for lifting
numEigvecs = 2;         % dimension of the embedding
weight = 5;             % weight for the diffusion map
foptions = optimoptions(@fsolve, 'Display', 'iter');% fsolve settings
options = odeset('AbsTol',10^-8,'RelTol',10^-8);    % ODE 45 settings

%% load diffusion map data
if full
    hways = readmatrix('../data/data_headways.csv');
    diffMap2D = DiffusionMap(hways, numEigvecs, weight);
    stepSize = .0025;       
    steps = 35;
else
    newData = readmatrix('../data/1000data2D.csv');
    diffMap2D = DiffusionMap(newData, numEigvecs, weight);
    stepSize = .01;      
    steps = 49;
end

%% load the reference states and previous bifurcation diagram
bif = readmatrix('../results/microBif.csv');
sigmas = std(bif(1:numCars,:));
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

%% load the reference states and comparison bifurcation diagram
start = 1;
change = 5;
v0_base2 = vel(start + change);
v0_base1 = vel(start);
T1 = T(start);
T2 = T(start + change);
p2 = norm(embed(:, start + change));
p1 = norm(embed(:, start));
coord2 = [p2; T2; v0_base2];
coord1 = [p1; T1; v0_base1];
sigma_1 = sigmas(start);
sigma_2 = sigmas(start + change);

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
bif = zeros(4,steps + 2);
bif(:, 1) = [coord1; sigma_1];
bif(:, 2) = [coord2; sigma_2];

%% pseudo arc length continuation
for iEq=3:steps+2
    fprintf('Starting iteration %d of %d \n', iEq-2, steps);
    
    % get the initial guess
    w = coord2 - coord1;
    scaledW = w ./ coord2;
    newGuess = coord2 + stepSize *(w/norm(w)); % first guess on the secant line
    scaledNewGuess = newGuess ./ coord2;
    
    % solve for the point on the curve
    u = fsolve(@(u)FW(u, diffMap2D, coord2, scaledW, scaledNewGuess), newGuess, foptions);
    [p, bif(4,iEq)] = ler( [u(1); 0], diffMap2D, u(3));
    bif(1, iEq) = norm(p);
    bif(2:3,iEq) = u(2:3)';
    
    %%reset the values for the arc length continuation
    coord1 = coord2;
    coord2 = u;
    
    % plot the diagram
    subplot(1,2,1); hold on;
    scatter(bif(3, iEq), bif(1,iEq), 400, 'b.'); drawnow;
    subplot(1,2,2); hold on;
    scatter(bif(1,iEq), bif(2,iEq), 400, 'b.'); drawnow;
end

% save results
if full
    writematrix(bif, '../results/bifurcation2D.csv');
else
    writematrix(bif, '../results/1000bifurcation2D.csv');
end


%% shooting part of the continuation
% PARAMETERS:
% u         - initial coordinate
% diffMap2D - diffusion map object
% coord2    - secant coordinate for normalization
% RETURNS:
% out  - the distance between the input point and the point determined by
%               lifting, evolving for T time, and restricting and dot
%               product with the secant line
    function out = periodicDistance(u, v, diffMap2D)
        
        inputPoint = [u(1);0];
        timeEvolve = nPeriod*u(2);
        [r0, ~, r1] = ler(inputPoint, diffMap2D, v, timeEvolve);
        
        % convert to polar coordinates
        angle0 = atan2(r0(2), r0(1));
        angle1 = atan2(r1(2), r1(1));
        out = [norm(r1) - norm(r0), angle1-angle0];
    end

%% function to zero
% PARAMETERS:
% u         - initial coordinate
% diffMap2D - diffusion map object
% coord2    - secant coordinate for normalization
% scaledW   - slope of the guess scaled by coord2
% scaledNewGuess    - prediction scaled by coord2
% RETURNS:
% F  - the distance between the input point and the point determined by
%               lifting, evolving for T time, and restricting and dot
%               product with the secant line
    function F = FW(u, diffMap2D, coord2, scaledW, scaledNewGuess)
        F = zeros(1,3);
        F(1:2) = periodicDistance(u(1:2), u(3), diffMap2D);
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
% coord      - the macrovariable after evolving for tskip
% sigma      - the std of the headways after evolving for tskip
% coord2      - the macrovariable  after evolving for tskip + nPeriod*T
    function [coord, sigma, coord2] = ler(newval, diffMap2D, v0, t)
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