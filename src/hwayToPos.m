% takes headways of cars and returns their positions
function positions = hwayToPos(headways)
numCars = 30;
positions = zeros(size(headways));
for ih = 2:numCars
    positions(ih) = positions(ih - 1) + headways(ih);
end
positions = circshift(positions, 1); % shift back
end