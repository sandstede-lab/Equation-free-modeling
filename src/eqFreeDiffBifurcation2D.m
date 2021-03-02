function eqFreeDiffBifurcation2D()
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 30;           % number of cars
steps = 100;                                % number of steps to take around the curve
stepSize = .01;                            % step size for the secant line approximation
full = true;
tskip = 500;
nPeriod = 1;
k = 8;
debugPlot = false;
foptions = optimoptions(@fsolve, 'TolFun',1e-15,'TolX',1e-15, 'Display', 'iter');

%% load diffusion map data
if full
    load('../data/diffMap2D.mat', 'diffMap2D');
else
    load('../data/1000diffMap2D.mat', 'diffMap2D');
end

% load the reference states and previous bifurcation diagram
load('../data/microBif.mat', 'bif');
vel = bif(end, :);
embed = diffMap2D.restrict(bif(1:numCars,:));
n = sqrt(embed(1,:).^2 + embed(2,:).^2);
start = 1;
v0_base2 = vel(start + 1);
v0_base1 = vel(start);
ref_2 = bif(1:numCars, start + 1);
ref_1 = bif(1:numCars, start);
T = -numCars./bif(numCars+1,:);
T1 = T(start);
T2 = T(start + 1);

embed_2 = diffMap2D.restrict(ref_2);
embed_1 = diffMap2D.restrict(ref_1);
p2 = norm(embed_2);
p1 = norm(embed_1);
coord2 = [p2, T1, v0_base2];
coord1 = [p1, T2, v0_base1];
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
scatter(n, T, 50, 'r.');
scatter(p1,T1, 400, 'k.'); drawnow;
scatter(p2, T2, 400, 'k.'); drawnow;
xlabel('\rho');
ylabel('T');

%% initialize secant continuation
bif = zeros(3,steps);
sigma = zeros(1,steps);

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
    sigma(iEq) = sig; % save the standard deviation of the headways to compare
    
    if full
        save('../data/newtonContinuation2D.mat', 'bif',  'sigma');
    else
        save('../data/1000newtonContinuation2D.mat', 'bif', 'sigma');
    end
    
    % plot the diagram
    subplot(1,2,1); hold on;
    scatter(bif(3, iEq), u(1), 400, 'b.'); drawnow;
    subplot(1,2,2); hold on;
    scatter(u(1), u(2), 400, 'b.'); drawnow;
    
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
        sigma = std(evoCars);
        
        if nargin > 3
            [~,evo2] = ode45(@microsystem,[0 t],evo(end,:), options,[v0 len h]);
            evoCars2 = getHeadways(evo2(end, 1:numCars)',len);
            sigma = std(evoCars2);
            coord2 = diffMap2D.restrict(evoCars2);
            
            if debugPlot
                clf;
                subplot(1,2,1); hold on;
                plot(1:30, evoCars);
                plot(1:30, evoCars2);
                subplot(1,2,2); hold on;
                scatter(diffMap2D.evecs(:,1), diffMap2D.evecs(:,2));
                scatter(coord(1), coord(2), 200, 'r.');
                scatter(coord2(1), coord2(2), 200, 'k.');
                drawnow;
            end
        end
    end

end
