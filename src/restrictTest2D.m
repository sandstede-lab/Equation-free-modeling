%% load diffusion map data
load('../data/diffMap2D.mat', 'hways', 'eps', 'evecs', 'evals');
allData = hways;

numRestrict = size(evecs,1);
restricted = zeros(numRestrict, 2);
restrictDiff = zeros(numRestrict, 1);
percentError = zeros(numRestrict,1);
origRestricted = zeros(numRestrict,2);
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
    origRestricted(i,:) = evecs(i, :);
    evecs2 = evecs;
    evecs2(i,:) = [];  
    restricted(i,:) = diffMapRestrict(allData(:,i),evals,evecs2,allData2,eps); % restrict this point
    restrictDiff(i) = norm((origRestricted(i,:) - restricted(i,:)));
    percentError(i) = norm((origRestricted(i,:) - restricted(i,:))./origRestricted(i,:));
end

disp(mean(percentError));

figure;
scatter(restricted(:,1), restricted(:,2), 200, restrictDiff, '.');
colorbar;
title('Restricted Points Colored by Distance from Original Embedding', 'fontsize', 12);
xlabel('\psi_1', 'fontsize', 12);
ylabel('\psi_2', 'fontsize', 12);
