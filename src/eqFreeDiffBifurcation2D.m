 function eqFreeDiffBifurcation2D()
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 30;           % number of cars
lsqOptions = optimset('Display','iter'); %lsqnonlin options
options2 = optimoptions('lsqnonlin');
options2.OptimalityTolerance = 1e-10;
options2.FunctionTolerance = 1e-10;
options2.Display = 'iter';

%% load diffusion map data
load('../data/diffMap2D.mat', 'hways', 'eps', 'evecs', 'evals');
allData=hways;
addKernelData(allData, eps); % save kernel data to avoid computing every time

% load the reference states and previous bifurcation diagram 
load('../data/microBif.mat', 'bif', 'vel', 'n');
start = 110;                                % location on the curve to start at
change = 1;
v0_base2 = vel(start + change);
v0_base1 = vel(start);
ref_2 = bif(1:numCars,start + change);
ref_1 = bif(1:numCars,start);

T2 = -numCars/bif(numCars + 1, start + change);
T1 = -numCars/bif(numCars + 1, start);
embed_2 = diffMapRestrict(ref_2,evals,evecs,allData,eps);
embed_1 = diffMapRestrict(ref_1,evals,evecs,allData,eps);
p2 = norm(embed_2);
p1 = norm(embed_1);
coord2 = [p2, T2, v0_base2];
coord1 = [p1, T1, v0_base1];
rayAngle = atan2(embed_2(2), embed_2(1));

% draw the bifurcation diagrams
figure;
subplot(1,2,1);hold on;
scatter(vel(100:end), n(100:end), 50, 'r.');
xlabel('v_0');
ylabel('\rho');
scatter(v0_base1, p1, 400,'k.'); drawnow;
scatter(v0_base2, p2, 400,'k.'); drawnow;

subplot(1,2,2); hold on;
scatter(n(100:end), -numCars./bif(numCars + 1, 100:end), 50, 'r.');
scatter(p1,T1, 400, 'k.'); drawnow;
scatter(p2, T2, 400, 'k.'); drawnow;
xlabel('\rho');
ylabel('T');

%% initialize secant continuation
steps = 500;                                % number of steps to take around the curve
bif = zeros(3,steps);                       % array to hold the bifurcation values
sigma = zeros(1,steps);
guesses = zeros(3,steps);
stepSize = .003;                            % step size for the secant line approximation

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = coord2 - coord1;
    scaledW = w ./ coord2;
    newGuess = coord2 + stepSize *(w/norm(w)); % first guess on the secant line
    guesses(:,iEq) = newGuess;
   
    scaledNewGuess = newGuess ./ coord2;

    %% alternate Newton's method using lsqnonlin
    [u, ~, ~, exitFlag] = lsqnonlin(@(u)periodicDistance(u,allData,evecs,evals,eps,coord2,rayAngle,scaledW,...
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
% allData
% evecs
% evals
% lereps
% coord2
% ray
% scaledW
% scaledNewGuess
% RETURNS:
% distance  - the distance between the input point and the point determined by 
%               lifting, evolving for T time, and restricting
% onplane   - whether this is on the plane orthogonal to w at the starting guess
    function out = periodicDistance(u,allData,evecs,evals,lereps, coord2, ray, scaledW, scaledNewGuess)
        inputPoint = [u(1)*cos(ray) ; u(1)*sin(ray)];
        [finalPoint, sig] = ler(inputPoint,allData,25*u(2),u(3),evecs,evals,lereps);
        
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
    function [coord, sigma, evo] = ler(newval,orig,t,v0,eigvecs,eigvals,lereps)
        options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
        lifted = diffMapLift(newval, eigvecs, eigvals, lereps,v0, orig, h);
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,[v0 len h]);
        evo = evo(end,:)';
        evoCars = getHeadways(evo(1:numCars),len);
        sigma = std(evoCars);
        coord = diffMapRestrict(evoCars,eigvals,eigvecs, orig, lereps);
    end

end