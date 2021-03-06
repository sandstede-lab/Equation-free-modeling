function eqFreeDiffBifurcation()
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;              % length of the ring road
numCars = 30;           % number of cars
alignTo = 10;
tskip = 300;
delta = 350; 
options = odeset('AbsTol',10^-8,'RelTol',10^-8);    % ODE 45 options
foptions = optimoptions(@fsolve,'Display','iter', 'TolFun',1e-10,'TolX',1e-10);
full = false;
k = 3; % number of points to use for lifting 
numEigvecs = 1;
weight = 5;

%% load diffusion map data
if full
    alignData = readmatrix('../data/alignData.csv');
    diffMap1D = DiffusionMap(alignData, numEigvecs, weight);
    numSteps = 200;
    stepSize = 0.001;
else
    newAlignData = readmatrix('../data/1000diffMap1D.csv');
    diffMap1D = DiffusionMap(newAlignData, numEigvecs, weight);
    numSteps = 32;     
    stepSize = .005;
end

% load the reference states and comparison bifurcation diagram
bif = readmatrix('../results/microBif.csv');
vel = bif(end, :);

start = length(vel)-1;
change = -1;
v0_base2 = vel(start+change);
v0_base1 = vel(start);
ref_2 = bif(1:numCars,start+change);
ref_1 = bif(1:numCars,start);
sigma_2 = std(ref_2);
sigma_1 = std(ref_1);
ref_2 = alignMax(ref_2, alignTo);
ref_1 = alignMax(ref_1, alignTo);

%initial psi values for secant line approximation
psi_1 = diffMap1D.restrict(ref_1);
psi_2 = diffMap1D.restrict(ref_2);

% draw the reference bifurcation diagrams
oldBif = bif;
figure; hold on;
scatter(vel, std(oldBif(1:numCars, :)), 200, 'r.');
scatter(v0_base1, sigma_1, 400,'k.'); drawnow;
scatter(v0_base2, sigma_2, 400,'k.'); drawnow;

%% initialize secant continuation                              
bif = zeros(3,numSteps);                       % array to hold the bifurcation values

%% pseudo arc length continuation
for iEq=1:numSteps
    fprintf('Starting iteration %d of %d \n', iEq, numSteps);
    w = [psi_2 - psi_1 ; v0_base2 - v0_base1];          % slope of the secant line
    newGuess = [psi_2; v0_base2] + stepSize *(w/norm(w)); % first guess on the secant line
     
    u = fsolve(@(u)FW(u, diffMap1D ,w, newGuess), newGuess, foptions);
    [bif(1,iEq), sig] = ler(u(1), diffMap1D, u(2));
    scatter(u(2), sig, 'b*'); drawnow;                          % plot the bifurcation diagram as it grows
    
    %% reset the values for the arc length continuation
    psi_1 = psi_2;
    v0_base1 = v0_base2;
    v0_base2 = u(2);
    psi_2 = u(1);                                             % find the new reference state
    bif(2,iEq) = v0_base2;                          % save the new solution
    bif(3, iEq) = sig;
    if full
        writematrix(bif, '../results/bifurcation1D.csv');
    else
        writematrix(bif, '../results/1000bifurcation1D.csv');
    end
end

%% function to zero for fsolve
% u         - the current value of (sigma, v0) that we're trying to find
%               with Newton's method
% ref       - the most recent reference state
% W         - the slope of the secant line for arc length continuation
% newGuess  - the first guess on the secant line for arc length
%               continuation
%
% RETURNS:
% fw    - the functions F and w evaluated at these parameters
    function [fw, J] = FW(u, diffMap1D, W, newGuess)
        fw = zeros(2,1);
        fw(1) = F(diffMap1D,u(1),u(2));
        fw(2) = W(1)*(u(1) -newGuess(1)) + W(2)*(u(2) - newGuess(2));
        if(nargout >1)
            J = jacobian(diffMap1D, u(1), u(2), W);
        end
    end

%% lift, evolve, restrict
% newval - the current macrovariable used to seed the lifting
% diffMap1D - diffusion map.
% v0    - the optimal velocity parameter for this state
% tReference - time to evolve
% RETURNS:
% psi     - new macrovariable
% sigma   - standard deviation of the evolved headways
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