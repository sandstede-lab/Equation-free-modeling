%% load diffusion map data
load('../data/diffMap2D.mat', 'hways', 'eps', 'evecs', 'evals', 'vel');
allData=hways;

numRestrict = size(evecs,1);
restricted = zeros(numRestrict, 2);
restrictDiff = zeros(numRestrict, 1);
percentError = zeros(numRestrict,1);

for i=1:numRestrict
    if mod(i, 500)==0
        disp(i);
    end
    lifted = diffMapLift(evecs(i, :)', evecs, evals, eps, vel(iRestrict), allData, h);
    restricted(i,:) = diffMapRestrict(getHeadways(lifted(1:30), 60),evals,evecs,allData,eps);
    if(norm(restricted(i,:)) > .025)
        disp(i);
    end
    restrictDiff(i) = norm(restricted(i,:) - evecs(i,:));
    percentError(i) = norm((evecs(i,:) - restricted(i,:))./evecs(i,:));
end

disp(mean(percentError));

figure;
scatter(restricted(:,1), restricted(:,2), 200, restrictDiff, '.');
cb = colorbar();
cb.Ruler.Scale = 'log';
cb.Ruler.MinorTick = 'on';
title('Lifted and Restricted Points Colored by Distance from Original Embedding', 'fontsize', 16);
xlabel('\psi_1', 'fontsize', 18);
ylabel('\psi_2', 'fontsize', 18);
