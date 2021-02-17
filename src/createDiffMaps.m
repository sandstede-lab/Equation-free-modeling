clear;

%% load traffic data
load('../data/30data.mat', 'hways', 'allData')

rng(18);                    % set random seed
numReduce = 1000;           % number of points to reduce the diffusion map to
numData = dims(hways, 2);   % number of points in the dataset
numCars = dims(hways, 1);   % number of cars
alignTo = 10;               % car to align data to
weight = 5;                 % median weight to choose epsilon 
linearFit = false;          % compute the linear fit or not
plotMaps = true;            % plot the diffusion maps or not
numEigvecs = 30;            % number of evectors to return for the linear fit

%% 1D diffusion map

% align data for 1D map
alignData = hways;
for i = 1:numData
    alignData(:,i) = alignMax(alignData(:,i), alignTo);
end

% take 1 eigenvector if not doing the linear fit, otherwise compute many
if ~linearFit
    numEigvecs = 1;
end
diffMap1D = DiffusionMap(alignData, numEigvecs, weight);

%calculate how unique each eigen direction is for a small sample
if linearFit
    r1 = zeros(numEigvecs, 1);
    idx = randsample(numData, 1000);
    r1(1) = 1;
    for j = 2:numEigvecs
        r1(j) = linearFit(diffMap1D.evecs(idx, :),j);
    end

    % reduce diffusion map to 1D
    diffMap1D.evecs = diffMap1D.evecs(:,1);
    diffMap1D.evals = diffMap1D.evals(1);
end
save('../data/diffMap1D.mat', 'diffMap1D');

% plot std vs diff map embedding
if plotMaps
    figure;
    scatter(diffMap1D.evecs, std(diffMap1D.data), 100,'.');
    xlabel('\psi_1', 'FontSize',14);
    ylabel('\sigma', 'FontSize',14);
    title('\sigma vs. \psi_1','FontSize',14);
end

% sort the diffusion map and only keep 1000 points evenly spaced
[~, idx] = sort(diffMap1D.evecs);
keep = floor(linspace(1, numData, numReduce));
% subset data
sortedData = diffMap1D.data(:, idx);
newAlignData = sortedData(:, keep);

% rerun diffMap with these points
diffMap1D = DiffusionMap(newAlignData, numEigvecs, weight);
save('../data/1000diffMap1D.mat', 'diffMap1D');

% plot the new diffusion map
if plotMaps
    figure;
    scatter(diffMap1D.evecs, std(diffMap1D.data),100,'.')
    xlabel('\psi_1', 'FontSize',14);
    ylabel('\sigma', 'FontSize',14);
    title('\sigma vs. \psi_1','FontSize',14);
end

%% 2D diffusion map

% keep first 2 dimensions unless checking the linear fit
if ~linearFit
    numEigvecs = 2;
end
diffMap2D = DiffusionMap(hways, numEigvecs, weight);

%calculate how unique each eigen direction is for a small sample
if linearFit
    r2 = zeros(numEigvecs, 1);
    idx = randsample(numData, 1000);
    r2(1) = 1;
    for j = 2:numEigvecs
        r2(j) = linearFit(diffMap2D.evecs(idx, :),j);
    end

    % plot the linear fit coefficients
    figure;
    hold on;
    scatter(1:numEigvecs, r1, 200, '.');
    scatter(1:numEigvecs, r2, 20, '*');
    legend('Aligned', 'Not-aligned');
    xlabel('\psi_k', 'FontSize',14);
    ylabel('r_k', 'FontSize',14);
    title('Local Linear Fit for Eigenvectors','FontSize',14);
    drawnow;

    % reduce diffusion map to 2D
    diffMap2D.evecs = diffMap2D.evecs(:,1:2);
    diffMap2D.evals = diffMap2D.evals(1:2, 1:2);
end
save('../data/diffMap2D.mat', 'diffMap2D');


[~, max1] = max(hways,[],1);  % locate the max headway for each data point
if plotMaps
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
    scatter(diffMap2D.evecs(:,1), diffMap2D.evecs(:,2), 100,  std(diffMap2D.data),'.');
    colorbar;
    color = colorbar;
    xlabel(color, '\sigma', 'fontsize', 14)
    colormap(jet);
    xlabel('\psi_1', 'FontSize',14);
    ylabel('\psi_2', 'FontSize',14);
    title('\psi_1 vs. \psi_2 Colored by Standard Deviation of the Headways','FontSize',14);
end

% convert to polar coordinates
n = sqrt(diffMap2D.evecs(:,1).^2 + diffMap2D.evecs(:,2).^2);
ray = atan2(diffMap2D.evecs(:,2), diffMap2D.evecs(:,1));

% sort the diffusion map by radial coordinate and only keep 2*nReduce points evenly spaced
[~,idx] = sort(n);
keep = floor(linspace(1, numData, 2*numReduce));
% subset angles
ray = ray(idx);
ray = ray(keep);
% subset data
sortedData = diffMap2D.data(:, idx);
newData = sortedData(:, keep);
% subset max
max1000 = max1(idx);
max1000 = max1000(keep);

% repeat with second coordinate (angular)
[~, idx] = sort(ray);
keep = floor(linspace(1, 2*numReduce, numReduce));
% subset data
sortedData = newData(:, idx);
newData = sortedData(:, keep);
% subset max
max1000 = max1000(idx);
max1000 = max1000(keep);

% rerun diff map with these points
diffMap2D = DiffusionMap(newData, numEigvecs, weight);
save('../data/1000diffMap2D.mat', 'diffMap2D');

if plotMaps
    % plot eigenvector 1 vs eigenvector 2, colored by max headway location
    figure;
    scatter(diffMap2D.evecs(:,1), diffMap2D.evecs(:,2), 100, max1000,'.'); hold on;
    color = colorbar;
    xlabel(color, 'Wave Position', 'fontsize', 14)
    colormap(jet);
    xlabel('\psi_1', 'FontSize',14);
    ylabel('\psi_2', 'FontSize',14);
    title('\psi_1 vs. \psi_2 Colored by Locations of the Max Headways','FontSize',14);
    drawnow;

    % plot eigenvector 1 vs eigenvector 2, colored by standard deviation
    figure;
    scatter(diffMap2D.evecs(:,1), diffMap2D.evecs(:,2), 100,  std(diffMap2D.data),'.');
    colorbar;
    color = colorbar;
    xlabel(color, '\sigma', 'fontsize', 14)
    colormap(jet);
    xlabel('\psi_1', 'FontSize',14);
    ylabel('\psi_2', 'FontSize',14);
    title('\psi_1 vs. \psi_2 Colored by Standard Deviation of the Headways','FontSize',14);

end
