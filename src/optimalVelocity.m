%% optimal velocity function given in paper
%    h         -  a parameter to represent the optimal velocity of the
%               	car
%    headway   - the distance between this car and the car ahead of it
%    v0        - ideal goal speed of each driver
%    returns:
%    v         - the optimal velocity of this car, who will speed or slow to try
%                   to meet it
function v = optimalVelocity(h, headway,v0)
v = v0 * (tanh(headway - h) + tanh(h));
end