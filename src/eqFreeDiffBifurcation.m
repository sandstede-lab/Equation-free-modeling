% eqFreeDiffBifurcation(true)
% computes the bifrucation diagram for the traffic system using aligned
% headways and a one-dimensional diffusion map
% 
% full  - if true uses the 5000 point diffusion map
%       - if false uses the 1000 point diffusion map
function eqFreeDiffBifurcation(full)
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 30;           % number of cars
alignTo = 10;           % position to align cars
tskip = 300;            % tiime to evolve to slow manifold
delta = 240;            % time to evolve for finite difference
k = 3;                  % number of points to use for lifting
numEigvecs = 1;         % dimension of the embedding
weight = 5;             % weight for the diffusion map
stepSize = .0025;       % continuation step size
options = odeset('AbsTol',10^-8,'RelTol',10^-8);    % ODE 45 settings
foptions = optimoptions(@fsolve,'Display','iter', ... % fsolve settings
    'TolFun',1e-12,'TolX',1e-12);

%% load diffusion map data
if full
    alignData = readmatrix('../data/alignData.csv');
    diffMap1D = DiffusionMap(alignData, numEigvecs, weight);
    numSteps = 42;
else
    newAlignData = readmatrix('../data/1000data1D.csv');
    diffMap1D = DiffusionMap(newAlignData, numEigvecs, weight);
    numSteps = 55;
end

%% load the reference states and comparison bifurcation diagram
bif = readmatrix('../results/microBif.csv');
vel = bif(end, :);
refs =bif(1:numCars,:);
for i=1:size(bif,2)
    refs(:, i) =  alignMax(refs(:,i), alignTo);
end
embed = diffMap1D.restrict(refs);

% save embedding of micro system
if full
    writematrix([embed', vel'], '../results/micro1D.csv');
else
    writematrix([embed', vel'], '../results/1000micro1D.csv');
end

%% get initiial points on the curve
start = 1;
change = 5;
v0_base2 = vel(start+change);
v0_base1 = vel(start);
psi_1 = embed(start);
psi_2 = embed(start + change);
coord1 = [psi_1; v0_base1];
coord2 = [psi_2; v0_base2];
sigma_1 = std(refs(:,start));
sigma_2 = std(refs(:,start + change));

% draw the reference bifurcation diagrams
figure; hold on;
scatter(vel, embed, 200, 'r.');
scatter(v0_base1, psi_1, 400,'k.'); drawnow;
scatter(v0_base2, psi_2, 400,'k.'); drawnow;

% initialize secant continuation
bif = zeros(3,numSteps + 2);
bif(:, 1) = [coord1; sigma_1];
bif(:, 2) = [coord2; sigma_2];

%% pseudo arc length continuation
for iEq=3:numSteps+2
    fprintf('Starting iteration %d of %d \n', iEq-2, numSteps);
    
    % get the initial guess
    w = coord2 - coord1;        
    wScaled = w ./ coord2;
    newGuess = coord2 + stepSize *(w/norm(w));
    scaledNewGuess = newGuess ./ coord2;
    
    % solve for the next point
    u = fsolve(@(u)FW(u, diffMap1D ,wScaled, scaledNewGuess, coord2), newGuess, foptions);
    bif(1:2, iEq) = u;
    [bif(1,iEq), bif(3,iEq)] = ler(u(1), diffMap1D, u(2));
    scatter(u(2), bif(1,iEq), 400, 'b.'); drawnow;
    
    % reset the values for the arc length continuation
    coord1 = coord2;
    coord2 = u;
end

%% save the results
if full
    writematrix(bif, '../results/bifurcation1D.csv');
else
    writematrix(bif, '../results/1000bifurcation1D.csv');
end

%% function to zero for fsolve
% u         - the current value of (sigma, v0) that we're trying to find
%               with Newton's method
% ref       - the most recent reference state
% W         - the slope of the secant line for arc length continuation
% newGuess  - the first guess on the secant line for arc length
%               continuation
% oldU      - previious u used to scale u for the secant dot product
% RETURNS:
% fw    - the functions F and w evaluated at these parameters
% J     - the jacobian evaluated at these parameters
    function [fw, J] = FW(u, diffMap1D, W, newGuess, oldU)
        fw = zeros(2,1);
        fw(1) = F(diffMap1D,u(1),u(2));
        scaledU = u ./oldU;
        fw(2) = dot(W, scaledU - newGuess);
        if(nargout >1)
            J = jacobian(diffMap1D, u(1), u(2), W);
        end
    end

%% lift, evolve, restrict
% newval    - the current macrovariable used to seed the lifting
% diffMap1D - diffusion map.
% v0        - the optimal velocity parameter for this state
% tReference - time to evolve
% RETURNS:
% psi     - new macrovariable after tskip
% sigma   - standard deviation of the evolved headways
% psi2    - new marcrovariable after tskip + delta
    function [psi, sigma, psi2] = ler(newval,diffMap1D,v0,tReference)
        liftedHeadway = diffMap1D.lift(newval, k);
        evoCars = [hwayToPos(liftedHeadway) ; optimalVelocity(h, liftedHeadway, v0)];
        
        [~,evo] = ode45(@microsystem,[0 tskip],evoCars,options,[v0 len h]);
        evoCars = getHeadways(evo(end, 1:numCars)', len);
        sigma = std(evoCars);
        psi = diffMap1D.restrict(alignMax(evoCars, alignTo));
        
        if nargin > 3
            [~,evo2] = ode45(@microsystem,[0 tReference],evo(end,:),options,[v0 len h]);
            evoCars2 = getHeadways(evo2(end, 1:numCars)', len);
            sigma = std(evoCars2);
            psi2 = diffMap1D.restrict(alignMax(evoCars2, alignTo));
        end
    end

%% Finite Difference Quotient
%  diffMap1D - diffusion map
%  sigma - the current value of the macrovariable
%  v0 - the velocity parameter for this state
%  RETURNS:
%  dif - the difference which approximates the time derivative
    function dif = F(diffMap1D, sigma, v0)
        [r0, ~, r1] = ler(sigma, diffMap1D, v0, delta);
        dif = (r1-r0)/delta;
    end

%% Jacobian for newton's method
% ref - The previous reference state used to compute F.
% sigma - the current value of sigma
% v0 - the velocity parameter for this state
% w - the secant direction
% J- The Jacobian, which will be given by
% | F_sigma    F_v0  |
% | w_sigma    w_vo  |
    function J = jacobian(diffMap1D, sigma, v0, w)
        J = zeros(2);
        delSigma = 0.00001;     % delta sigma used for finite difference of F
        delv0 = 0.00001;        % delta v0 used for finite difference of F
        unchanged = F(diffMap1D, sigma, v0);
        J(1,1) = (F(diffMap1D, sigma + delSigma, v0) - unchanged)/delSigma;
        J(1,2) = (F(diffMap1D, sigma,v0 + delv0) - unchanged)/delv0;
        J(2,:) = w';
    end

end