h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
alignTo = 20;
tskip = 300;            % times for evolving
delta = 2000;

stepSize = .0001;        % step size for the secant line approximation
delSigma = 10^(-7);     % delta sigma used for finite difference of F
delv0 = 10^(-7);        % delta v0 used for finite difference of F
tolerance = 10^(-14);    % tolerance for Newton's method

options = odeset('AbsTol',10^-8,'RelTol',10^-8);                                % ODE 45 options
foptions = optimoptions(@fsolve, 'TolFun', tolerance, 'display', 'off');      	%fsolve options

%% load diffusion map data
load('moreData.mat', 'betterData');
allData = betterData;                                   % aligned headway data

numEigvecs = 1;                                         % number of eigenvectors to return
[evecs, evals, eps] = runDiffMap(allData,numEigvecs);   % run the diffusion map

[evecs,ia,~] = unique(evecs);
allData = allData(:, ia);

%{
% calculate how unique each eigen direction is
r = zeros(numEigvecs, 1);
r(1) = 1;
for j = 2:numEigvecs
    disp(j);
    r(j) = linearFit(evecs,j);
end
%}

%% initialize secant continuation
steps = 400;                                % number of steps to take around the curve
bif = zeros(2,steps);                       % array to hold the bifurcation values

% initialize the first reference state
load('start0.884800.mat', 'trafficData', 'vel');
ref_2 = getHeadways(trafficData(1:numCars), len);
ref_2 = alignMax(ref_2, alignTo, false);
v0_base2 = vel;

% initialize the second reference state
load('start0.885000.mat', 'trafficData', 'vel');
ref_1 = getHeadways(trafficData(1:numCars), len);
ref_1 = alignMax(ref_1, alignTo, false);
v0_base1 = vel;

sigma_1 = diffMapRestrict(ref_1,evals,evecs,allData,eps);      %initial sigma values for secant line approximation
sigma_2 = diffMapRestrict(ref_2,evals,evecs,allData,eps);

%% pseudo arc length continuation
figure;
hold on;
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = [sigma_2 - sigma_1 ; v0_base2 - v0_base1];          % slope of the secant line
    newGuess = [sigma_2; v0_base2] + stepSize *(w/norm(w)); % first guess on the secant line

    %% Newton's method using fsolve
    u = fsolve(@(u)FW(u,allData,w,newGuess,evecs,evals,eps), newGuess,foptions);
    fprintf('\t Velocity is: %f \n', u(2));
    sig = interp1(evecs,std(allData(1:numCars,:)),u(1));        % interpolate to sigma
    scatter(u(2), sig, 'b*'); drawnow;                          % plot the bifurcation diagram as it grows
    
    %% reset the values for the arc length continuation
    sigma_1 = sigma_2;
    v0_base1 = v0_base2;
    v0_base2 = u(2);
    sigma_2 = u(1);                                             % find the new reference state
    bif(:,iEq) = [sigma_2 ; v0_base2];                          % save the new solution
end
hold off;

% plot sigma vs eigenvector 1
figure;
scatter(std(allData(1:numCars, :)), evecs(:,1),'b.');
xlabel('\sigma');
ylabel('Steve');

%% plot the bifurcation diagram
figure;
scatter(bif(2,:),bif(1,:),'*');
xlabel('v0');
ylabel('\Phi_1');

% interpolate back to the standard deviation values
sig = interp1(evecs,std(allData(1:numCars,:)),bif(1,:));
% plot the bifurcation diagram using standard deviation coordinates
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
    function [fw, J] = FW(u,ref,W,newGuess,evecs,evals,lereps)
        fw = zeros(2,1);
        fw(1) = F(ref,u(1),u(2),evecs,evals,lereps);
        fw(2) = W(1)*(u(1)-newGuess(1)) + W(2)*(u(2) - newGuess(2));
        if(nargout >1)
            J = jacobian(ref, u(1), u(2),W,evecs,evals,lereps);
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
    function [sigma,new_state, sigma2, new_state2] = ler(newval,orig,t,v0,eigvecs,eigvals,lereps,tReference)
        lifted = smartLift(newval, eigvecs, eigvals, lereps,v0, orig, alignTo, h, len);
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,[v0 len h]);
        if (nargin > 7)
            [~,evo2] = ode45(@microsystem,[0 tReference],evo(end,1:2*numCars)',options,[v0 len h]);
            evo2Cars = evo2(end, :)';
            evo2Cars = getHeadways(evo2Cars(1:numCars), len); 
            evo2Cars = alignMax(evo2Cars, alignTo, false);
            sigma2 = diffMapRestrict(evo2Cars,eigvals,eigvecs, orig, lereps);
        end
        evoCars = evo(end, :)';
        evoCars = getHeadways(evoCars(1:numCars), len);
        evoCars = alignMax(evoCars, alignTo, false);
        sigma = diffMapRestrict(evoCars, eigvals, eigvecs, orig, lereps);
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
    function dif = F(ref, sigma,v0,eigvecs,eigvals,lereps)
        [r0, ~, r1] = ler(sigma, ref, tskip, v0, eigvecs,eigvals,lereps, delta);
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
    function J = jacobian(ref, sigma, v0,w,eigvecs,eigvals,lereps)
        J = zeros(2);
        unchanged = F(ref, sigma, v0,eigvecs,eigvals,lereps);
        J(1,1) = (F(ref, sigma + delSigma, v0,eigvecs,eigvals,lereps) - unchanged)/delSigma;
        J(1,2) = (F(ref, sigma,v0 + delv0,eigvecs,eigvals,lereps) - unchanged)/delv0;
        J(2,:) = w';
    end
