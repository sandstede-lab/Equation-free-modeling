function eqFreeDiffBifurcation()
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;              % length of the ring road
numCars = 30;           % number of cars
alignTo = 10;
tskip = 100;            % times for evolving
delta = 1000;
stepSize = 0.001;%.00001;        % step size for the secant line approximation
delSigma = 10^(-7);     % delta sigma used for finite difference of F
delv0 = 10^(-7);        % delta v0 used for finite difference of F
options = odeset('AbsTol',10^-8,'RelTol',10^-8);    % ODE 45 options
options2 = optimoptions('lsqnonlin');

%% load diffusion map data
load('../data/diffMap1D.mat', 'diffMap1D');

% load the reference states
load('../data/microBif.mat', 'bif');
vel = bif(end, :);
start = 50;                                % location on the curve to start at
change = 1;
v0_base2 = vel(start + change);
v0_base1 = vel(start);
ref_2 = bif(1:numCars,start + change);
ref_1 = bif(1:numCars,start);
ref_2 = alignMax(ref_2, alignTo);
ref_1 = alignMax(ref_1, alignTo);

%initial psi values for secant line approximation
psi_1 = diffMap1D.restrict(ref_1);      
psi_2 = diffMap1D.restrict(ref_2);

% draw the reference bifurcation diagrams
oldBif = bif;
figure; hold on;
scatter(vel(start:end), std(oldBif(1:numCars, start:end)), 200, 'r.');
sig = interp1(diffMap1D.evecs, std(diffMap1D.data(1:numCars,:)),psi_1);        % interpolate to sigma
scatter(v0_base1, sig, 400,'k.'); drawnow;
sig = interp1(diffMap1D.evecs, std(diffMap1D.data(1:numCars,:)),psi_2);        % interpolate to sigma
scatter(v0_base2, sig, 400,'k.'); drawnow;

%% initialize secant continuation
steps = 400;                                % number of steps to take around the curve
bif = zeros(2,steps);                       % array to hold the bifurcation values

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = [psi_2 - psi_1 ; v0_base2 - v0_base1];          % slope of the secant line
    newGuess = [psi_2; v0_base2] + stepSize *(w/norm(w)); % first guess on the secant line

    %% Newton's method using lsqnonlin
    u = lsqnonlin(@(u)FW(u, diffMap1D ,w, newGuess), newGuess,[min(diffMap1D.evecs) 0.9],[max(diffMap1D.evecs) 1.2],options2);
    sig = interp1(diffMap1D.evecs,std(diffMap1D.data(1:numCars,:)),u(1));        % interpolate to sigma
    scatter(u(2), sig, 'b*'); drawnow;                          % plot the bifurcation diagram as it grows
    
    %% reset the values for the arc length continuation
    psi_1 = psi_2;
    v0_base1 = v0_base2;
    v0_base2 = u(2);
    psi_2 = u(1);                                             % find the new reference state
    bif(:,iEq) = [psi_2 ; v0_base2];                          % save the new solution
    save('../data/newtonContinuation1D.mat', 'bif');
end

%% plot the bifurcation diagram
figure;
scatter(bif(2,:),bif(1,:),'*');
xlabel('v0');
ylabel('\Phi_1');

% interpolate back to the standard deviation values and plot
sig = interp1(diffMap1D.evecs,std(diffMap1D.data(1:numCars,:)),bif(1,:));
figure;
scatter(bif(2,:),sig);
title('Interpolated sigmas');
xlabel('v0');
ylabel('\sigma');

%% functions

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
    function [fw, J] = FW(u, diffMap1D, W,newGuess)
        fw = zeros(2,1);
        fw(1) = F(diffMap1D,u(1),u(2));
        fw(2) = W(1)*(u(1)-newGuess(1)) + W(2)*(u(2) - newGuess(2));
        if(nargout >1)
            J = jacobian(diffMap1D, u(1), u(2),W);
        end
    end

%% lift, evolve, restrict
% sigma - the current value of the std, used to seed the lifting
% ref   - a reference state, used to seed the lifting
% t     - the duration to evaluate the lifted parameters
% p     -  parameter which will degrade the accuracy of the lifting
%               operator.  Should be kept at p = 1 by default.
% v0    - the optimal velocity parameter for this state
% RETURNS:
% sigma     - the std. of the headways after restricting has occured
% new_state - the final state of the evolution, which can be used as a
%               future reference state
    function [sigma,new_state, sigma2, new_state2] = ler(newval,diffMap1D,t,v0,tReference)
        liftedHeadway = diffMap1D.lift(newval);
        lifted = [hwayToPos(liftedHeadway) ; optimalVelocity(h, liftedHeadway, v0)];
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,[v0 len h]);
        if (nargin > 4)
            [~,evo2] = ode45(@microsystem,[0 tReference],evo(end,1:2*numCars)',options,[v0 len h]);
            evo2Cars = evo2(end, :)';
            evo2Cars = getHeadways(evo2Cars(1:numCars), len); 
            evo2Cars = alignMax(evo2Cars, alignTo);
            sigma2 = diffMap1D.restrict(evo2Cars);
        end
        evoCars = evo(end, :)';
        evoCars = getHeadways(evoCars(1:numCars), len);
        evoCars = alignMax(evoCars, alignTo);
        sigma = diffMap1D.restrict(evoCars);
        if(nargout >= 2)
            new_state = evo(end,1:2*numCars)';
        end
        if(nargout == 4)
            new_state2 = evo2(end,1:2*numCars)';
        end
    end

%% Finite Difference Quotient
%  ref - reference state to base the lifting
%  sigma - the current value of the std. at which point to approximate
%       the time derivative
%  v0 - the velocity parameter for this state
%  RETURNS:
%  dif - the difference which approximates the time derivative
    function dif = F(diffMap1D, sigma, v0)
        [r0, ~, r1] = ler(sigma, diffMap1D, tskip, v0, delta);
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
        unchanged = F(diffMap1D, sigma, v0);
        J(1,1) = (F(diffMap1D, sigma + delSigma, v0) - unchanged)/delSigma;
        J(1,2) = (F(diffMap1D, sigma,v0 + delv0) - unchanged)/delv0;
        J(2,:) = w';
    end

end