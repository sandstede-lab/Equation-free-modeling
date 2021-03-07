clear;

%% load traffic data
allData = readmatrix('../data/data_full.csv');
hways = readmatrix('../data/data_headways.csv');

figure;
scatter(allData(:,end), std(hways));

rng(18);                    % set random seed
numReduce = 1000;           % number of points to reduce the diffusion map to
numData = size(hways, 2);   % number of points in the dataset
numCars = size(hways, 1);   % number of cars
len = 60;                   % length of the ring road
alignTo = 10;               % car to align data to
weight = 5;                 % median weight to choose epsilon 
doLinearFit = true;        % compute the linear fit or not
plotMaps = true;           % plot the diffusion maps or not
numEigvecs = 30;            % number of evectors to return for the linear fit

%% 1D diffusion map

% align data for 1D map
alignData = hways;
for i = 1:numData
    alignData(:,i) = alignMax(alignData(:,i), alignTo);
end
writematrix(alignData, '../data/alignData.csv');

% take 1 eigenvector if not doing the linear fit, otherwise compute many
if ~doLinearFit
    numEigvecs = 1;
end
diffMap1D = DiffusionMap(alignData, numEigvecs, weight);

%calculate how unique each eigen direction is for a small sample
if doLinearFit
    r1 = zeros(numEigvecs, 1);
    idx = randsample(numData, 1000);
    r1(1) = 1;
    for j = 2:numEigvecs
        r1(j) = linearFit(diffMap1D.evecs(idx, :),j);
    end

    % plot the linear fit coefficients
    figure;
    hold on;
    scatter(1:numEigvecs, r1, 200, 'b.');
    xlabel('\psi_j', 'FontSize',14);
    ylabel('r_j', 'FontSize',14);
    title('Local Linear Fit for Eigenvectors for Full 1D Diffusion Map','FontSize',14);
    drawnow;
    
    % save r_j values
    writematrix(r1, '../results/r1.csv');

    % reduce diffusion map to 1D
    diffMap1D.evecs = diffMap1D.evecs(:,1);
    diffMap1D.evals = diffMap1D.evals(1);
end
writematrix(diffMap1D.evecs, '../results/embedding1D.csv');

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

figure;
vel = allData(idx, end);
vel = vel(keep);
scatter(vel, std(newAlignData));


%% rerun diffMap with these points
if doLinearFit
    numEigvecs = 30;
end
diffMap1D = DiffusionMap(newAlignData, numEigvecs, weight);

%calculate how unique each eigen direction is
if doLinearFit
    r1 = zeros(numEigvecs, 1);
    r1(1) = 1;
    for j = 2:numEigvecs
        r1(j) = linearFit(diffMap1D.evecs,j);
    end

    % plot the linear fit coefficients
    figure;
    hold on;
    scatter(1:numEigvecs, r1, 200, 'b.');
    xlabel('\psi_j', 'FontSize',14);
    ylabel('r_j', 'FontSize',14);
    title('Local Linear Fit for Eigenvectors for Reduced 1D Diffusion Map','FontSize',14);
    drawnow;
    
    % save r_j values
    writematrix(r1, '../results/r1_reduced.csv');

    % reduce diffusion map to 1D
    diffMap1D.evecs = diffMap1D.evecs(:,1);
    diffMap1D.evals = diffMap1D.evals(1);
end
writematrix(newAlignData, '../data/1000diffMap1D.csv');
writematrix(diffMap1D.evecs, '../results/1000embedding1D.csv');

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
if ~doLinearFit
    numEigvecs = 2;
end
diffMap2D = DiffusionMap(hways, numEigvecs, weight);

%calculate how unique each eigen direction is for a small sample
if doLinearFit
    r2 = zeros(numEigvecs, 1);
    idx = randsample(numData, 1000);
    r2(1) = 1;
    for j = 2:numEigvecs
        r2(j) = linearFit(diffMap2D.evecs(idx, :),j);
    end

    % plot the linear fit coefficients
    figure;
    hold on;
    scatter(1:numEigvecs, r2, 200, 'b.');
    xlabel('\psi_j', 'FontSize',14);
    ylabel('r_j', 'FontSize',14);
    title('Local Linear Fit for Eigenvectors for Full 2D Diffusion Map','FontSize',14);
    drawnow;
    
    % save r_j values
    writematrix(r2, '../results/r2.csv');

    % reduce diffusion map to 2D
    diffMap2D.evecs = diffMap2D.evecs(:, 1:2);
    diffMap2D.evals = diffMap2D.evals(1:2,1:2);
end
writematrix(diffMap2D.evecs, '../results/embedding2D.csv');

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
    scatter(diffMap2D.evecs(:,1), diffMap2D.evecs(:,2), 100,  std(hways),'.');
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

% sort the diffusion map by radial coordinate and only keep some points evenly spaced
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

%% rerun diff map with these points
if doLinearFit
    numEigvecs = 30;
end
diffMap2D = DiffusionMap(newData, numEigvecs, weight);

% compute local linear fit for reduced map
if doLinearFit
    r2 = zeros(numEigvecs, 1);
    r2(1) = 1;
    for j = 2:numEigvecs
        r2(j) = linearFit(diffMap2D.evecs,j);
    end

    % plot the linear fit coefficients
    figure;
    hold on;
    scatter(1:numEigvecs, r2, 200, 'b.');
    xlabel('\psi_j', 'FontSize',14);
    ylabel('r_j', 'FontSize',14);
    title('Local Linear Fit for Eigenvectors for Reduced 2D Diffusion Map','FontSize',14);
    drawnow;

    % save r_j values
    writematrix(r2, '../results/r2_reduced.csv');
    
    % reduce diffusion map to 2D
    diffMap2D.evecs = diffMap2D.evecs(:, 1:2);
    diffMap2D.evals = diffMap2D.evals(1:2, 1:2);
end
writematrix(newData, '../data/1000diffMap2D.csv');
writematrix(diffMap2D.evecs, '../results/1000embedding2D.csv');

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


