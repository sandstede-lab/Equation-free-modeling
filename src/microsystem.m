%% ODE that governs individual cars
% Governs the movement of the individual cars (microvariables)
% ~         - dummy parameter for time, to allow use in ode45
% params    -  v0 (optimal veloctiy parameter), len (length of the ring
%               road), and h (optimal velocity parameter)
% Returns:
% u         - a column vector of the distribution of the cars and their
%               velocities, of size 2*numCars.  The position of car i
%               and its velocity are given at num params(i),
%               params(i + numcars)
%
function u = microsystem(~,colCars, param)
v0 = param(1);
len = param(2);
h = param(3);
numCars = length(colCars)/2;
invT = 1.7;
headways = getHeadways(colCars(1:numCars), len);

u = zeros(2*numCars,1);
u(1:numCars,1) = colCars(numCars+1:2*numCars,1);
u(numCars+1:2*numCars,1) = invT*(optimalVelocity(h, headways,v0) - colCars(numCars+1:2*numCars,1));
end
