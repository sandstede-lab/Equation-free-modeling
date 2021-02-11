%% tests how accurate the restriction operator is when a data point is
% removed from the diffusion map

%% load diffusion map data
h=2.4;
load('../data/diffMap1D.mat', 'alignData', 'evals', 'evecs', 'eps', 'vel'); 
allData = alignData;

numRestrict = length(evecs);
orig = zeros(numRestrict, 1);                   % to store the data
restricted = zeros(numRestrict, 1);
restrictDiff = zeros(numRestrict, 1);
origRestricted = zeros(numRestrict,1);

percentError = zeros(numRestrict,1);
%nystrom error comparison, calculated the following way:
%   'correct' coordinates are given by the evec coord
%   approximated coordinates are given by diff map restrict, when the value
%   is removed from the data set
%   better approximated coordinates are given by diff map restrict, when
%   the value is still included in the data set
for i = 1:numRestrict
    if mod(i, 500)==0
        disp(i);
    end
    allData2 = allData;
    allData2(:, i) = [];                    % remove this point from the data set
    origRestricted(i) = evecs(i, :);
    evecs2 = evecs;
    evecs2(i) = []; 
    restricted(i) = diffMapRestrict(allData(:,i),evals,evecs2,allData2,eps); % restrict this point
    restrictDiff(i) = norm((origRestricted(i) - restricted(i)));
    percentError(i) = restrictDiff(i)/norm(origRestricted(i));
end

disp(mean(percentError));

figure; % plot the differences
scatter(evecs,restricted, 200, restrictDiff, '.');
colorbar;
title('Restricted Points Colored by Distance from Original Embedding', 'fontsize', 12);
xlabel('Original Coordinate', 'fontsize', 12);
ylabel('Restricted Coordinate', 'fontsize', 12);
