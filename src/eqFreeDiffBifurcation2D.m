 function eqFreeDiffBifurcation2D()
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 30;           % number of cars
steps = 800;                                % number of steps to take around the curve
stepSize = .003;                            % step size for the secant line approximation

lsqOptions = optimset('Display','iter'); %lsqnonlin options
options2 = optimoptions('lsqnonlin');
options2.OptimalityTolerance = 1e-10;
options2.FunctionTolerance = 1e-10;
options2.Display = 'iter';

%% load diffusion map data
load('../data/diffMap2D.mat', 'diffMap2D');

% load the reference states and previous bifurcation diagram 
load('../data/microBif.mat', 'bif');
vel = bif(end, :);
n = diffMap2D.restrict(bif(1:numCars,:));
n = sqrt(n(1,:).^2 + n(2,:).^2);
v0_base2 = vel(2);
v0_base1 = vel(1);
ref_2 = bif(1:numCars, 2);
ref_1 = bif(1:numCars, 1);
T2 = -numCars/bif(numCars + 1, 2);
T1 = -numCars/bif(numCars + 1, 1);

embed_2 = diffMap2D.restrict(ref_2);
embed_1 = diffMap2D.restrict(ref_1);
p2 = norm(embed_2);
p1 = norm(embed_1);
coord2 = [p2, T2, v0_base2];
coord1 = [p1, T1, v0_base1];
rayAngle = atan2(embed_2(2), embed_2(1));

% draw the bifurcation diagrams
figure;
subplot(1,2,1);hold on;
scatter(vel, n, 50, 'r.');
xlabel('v_0');
ylabel('\rho');
scatter(v0_base1, p1, 400,'k.'); drawnow;
scatter(v0_base2, p2, 400,'k.'); drawnow;

subplot(1,2,2); hold on;
scatter(n, -numCars./bif(numCars + 1, :), 50, 'r.');
scatter(p1,T1, 400, 'k.'); drawnow;
scatter(p2, T2, 400, 'k.'); drawnow;
xlabel('\rho');
ylabel('T');

%% initialize secant continuation
bif = zeros(3,steps);                       % array to hold the bifurcation values
sigma = zeros(1,steps);
guesses = zeros(3,steps);

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = coord2 - coord1;
    scaledW = w ./ coord2;
    newGuess = coord2 + stepSize *(w/norm(w)); % first guess on the secant line
    guesses(:,iEq) = newGuess;
    scaledNewGuess = newGuess ./ coord2;

    %% alternate Newton's method using lsqnonlin
    [u, ~, ~, exitFlag] = lsqnonlin(@(u)periodicDistance(u,diffMap2D,coord2,rayAngle,scaledW,...
        scaledNewGuess), newGuess,[],[],lsqOptions);
    
    bif(:,iEq) = u';                                            % save the new solution
    
    %% reset the values for the arc length continuation
    coord1 = coord2;
    coord2 = u;
    
    sigma(iEq) = sig; % save the standard deviation of the headways to compare
    
    fprintf('Exit flag: %d \n', exitFlag);
    save('../data/newtonContinuation2D.mat', 'bif', 'guesses', 'rayAngle', 'sigma');
    
    subplot(1,2,1); hold on;
    scatter(u(3), u(1), 400, 'b.'); drawnow;
    subplot(1,2,2); hold on;
    scatter(u(1), u(2), 400, 'b.'); drawnow;
end


%% function to zero
% PARAMETERS:
% u         - the pair of radius and period (sigma,T,v0) that we're trying to find with
% diffMap2D
% coord2
% ray
% scaledW
% scaledNewGuess
% RETURNS:
% distance  - the distance between the input point and the point determined by 
%               lifting, evolving for T time, and restricting
% onplane   - whether this is on the plane orthogonal to w at the starting guess
    function out = periodicDistance(u,diffMap2D, coord2, ray, scaledW, scaledNewGuess)
        inputPoint = [u(1)*cos(ray) ; u(1)*sin(ray)];
        [finalPoint, sig] = ler(inputPoint, diffMap2D, 25*u(2),u(3));
        
        angle = atan2(finalPoint(2), finalPoint(1));
        r = norm(finalPoint);
        distance = [((u(1)-r)/10^(-3))^2 ; (angle-ray)^2 ]; % compute the angular and radial distances

        scaledU = u ./ coord2;
        onPlane = dot(scaledW, scaledU - scaledNewGuess); % angle from the secant line
        out = [distance' onPlane];
    end


%% lift, evolve, restrict
% PARAMETERS:
% newval     - the current value of the macrovariables,
%                 used to seed the lifting
% orig       - the data set, for use in lifting
% t          - the duration to evaluate the lifted parameters
% v0         - the optimal velocity parameter for this state
% eigvecs    - eigenvectors from diffusion map
% eigvals    - eigenvalues from diffusion map
% lereps     - epsilon used in diffusionmap
% RETURNS:
% sigma      - the macro-state of the headways after evolving for t
    function [coord, sigma, evo] = ler(newval, diffMap2D, t,v0)
        options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
        liftedHeadway = diffMap2D.lift(newval);
        lifted = [hwayToPos(liftedHeadway) ; optimalVelocity(h, liftedHeadway, v0)];
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,[v0 len h]);
        evo = evo(end,:)';
        evoCars = getHeadways(evo(1:numCars),len);
        sigma = std(evoCars);
        coord = diffMap2D.restrict(evoCars);
    end

end
