%% tests the identity R(L) = I by taking coordinates in the diffusion map 
% embedding, lifting them, restricting them, and looking at the difference

% load diffusion map data
h=2.4;
load('../data/diffMap1D.mat', 'alignData', 'evals', 'evecs', 'eps', 'vel'); 

numRestrict = length(evecs);
orig = zeros(numRestrict, 1);
restricted = zeros(numRestrict, 1);
restrictDiff = zeros(numRestrict, 1);
percentError = zeros(numRestrict,1);

for i=1:numRestrict
    if mod(i, 500)==0
        disp(i);
    end
    orig(i) = evecs(i);
    lifted = diffMapLift(orig(i), evecs, evals, eps, vel(i), alignData, h);  % lift the profile
    restricted(i) = diffMapRestrict(getHeadways(lifted(1:30), 60),evals,evecs,alignData,eps); % restrict the lifted profile
    restrictDiff(i) = norm(restricted(i) - orig(i));    % compute the difference from the original embedding
    percentError(i) = abs(restricted(i)-orig(i))/abs(orig(i));
end

disp(mean(percentError));

figure; % plot the differences
scatter(orig,restricted, 200, restrictDiff, '.');
cb = colorbar();
cb.Ruler.Scale = 'log';
cb.Ruler.MinorTick = 'on';
title('Lifted and Restricted Points Colored by Distance from Original Embedding', 'fontsize', 16);
xlabel('Original Coordinates', 'fontsize', 18);
ylabel('New Coordinates', 'fontsize', 18);