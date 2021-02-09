%% function to calculate headways of car vector
%  v   - column vector of the cars' positions
%  len - the length of the ring road
%  Returns:
%  hways - column vector of cars', len headways
function hways = getHeadways(v, len)
futureCars = circshift(v,[-1,0]);
hways = mod(futureCars - v, len);
end