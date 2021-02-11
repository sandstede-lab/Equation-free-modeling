%% load traffic data
load('../data/30data.mat', 'hways', 'allData')
vel = allData(:,end);
dims = size(hways);
numCars = dims(1);
numData = dims(2);
alignTo = 10;

% align data for 1D map
alignData = hways;
for i = 1:numData
    alignData(:,i) = alignMax(alignData(:,i), alignTo);
end

%% 1D diffusion map

numEigvecs = 1;       % number of eigenvectors to return
weight = 5;  % median weight to choose epsilon 
diffMap1D = DiffusionMap(alignData, numEigvecs, weight);

save('../data/diffMap1D.mat', 'diffMap1D',  'vel');

figure;
scatter(diffMap1D.evecs, std(hways),100,'.');
xlabel('\psi_1', 'FontSize',14);
ylabel('\sigma', 'FontSize',14);
title('\sigma vs. \psi_1','FontSize',14);

%% 2D diffusion map

numEigvecs = 2;       % number of eigenvectors to return
diffMap2D = DiffusionMap(hways, numEigvecs, weight);

save('../data/diffMap2D.mat', 'diffMap2D', 'vel');

%{
%calculate how unique each eigen direction is
r = zeros(numEigvecs, 1);
r(1) = 1;
for j = 2:numEigvecs
    r(j) = linearFit(evecs,j);
end
%}

[~, max1] = max(hways,[],1);  % locate the max headway for each data point
% plot eigenvector 1 vs eigenvector 2, colored by max headway location
figure;
scatter(diffMap2D.evecs(:,1), diffMap2D.evecs(:,2), 100, max1,'.'); hold on;
color = colorbar;
xlabel(color, 'Wave Position', 'fontsize', 14)
colormap(jet);
xlabel('\psi_1', 'FontSize',14);
ylabel('\psi_2', 'FontSize',14);
title('\psi_1 vs. \psi_2 Colored by Locations of the Max Headways','FontSize',14);
drawnow;

% plot eigenvector 1 vs eigenvector 2, colored by standard deviation
figure;
scatter(diffMap2D.evecs(:,1), diffMap2D.evecs(:,2), 100,  std(hways),'.');
colorbar;
color = colorbar;
xlabel(color, '\sigma', 'fontsize', 14)
colormap(jet);
xlabel('\psi_1', 'FontSize',14);
ylabel('\psi_2', 'FontSize',14);
title('\psi_1 vs. \psi_2 Colored by Standard Deviation of the Headways','FontSize',14);
