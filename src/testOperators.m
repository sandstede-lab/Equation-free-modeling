%% tests how accurate the lifting/restricting operators are

%% load 1D diffusion map data
alignData = readmatrix('../data/alignData.csv');
diffMap1D = DiffusionMap(alignData, 1, 5);

% test the restriction operator
[percentError, restricted, restrictDiff] = diffMap1D.testRestrict();
disp(mean(percentError)); 

figure; % plot the differences
scatter(diffMap1D.evecs, restricted, 200, restrictDiff, '.');
colorbar;
title('Restricted Points Colored by Distance from Original Embedding', 'fontsize', 12);
xlabel('Original Coordinate', 'fontsize', 12);
ylabel('Restricted Coordinate', 'fontsize', 12);

% test the lift operator 
[percentError, restricted, restrictDiff] = diffMap1D.testLift(3);
disp(mean(percentError)); 

figure; % plot the differences
scatter(diffMap1D.evecs, restricted, 200, restrictDiff, '.');
cb = colorbar();
cb.Ruler.Scale = 'log';
cb.Ruler.MinorTick = 'on';
title('Lifted and Restricted Points Colored by Distance from Original Embedding', 'fontsize', 16);
xlabel('Original Coordinates', 'fontsize', 18);
ylabel('New Coordinates', 'fontsize', 18);

%% load 2D diffusion map data
 hways = readmatrix('../data/data_headways.csv');
 diffMap2D = DiffusionMap(hways, 2, 5);

% test the restriction operator
[percentError, restricted, restrictDiff] = diffMap2D.testRestrict();
disp(mean(percentError)); 

figure;
scatter(restricted(:,1), restricted(:,2), 200, restrictDiff, '.');
colorbar;
title('Restricted Points Colored by Distance from Original Embedding', 'fontsize', 12);
xlabel('\psi_1', 'fontsize', 12);

% test the lifting operator
[percentError, restricted, restrictDiff] = diffMap2D.testLift(8);
disp(mean(percentError)); 

figure;
scatter(restricted(:,1), restricted(:,2), 200, restrictDiff, '.');
cb = colorbar();
cb.Ruler.Scale = 'log';
cb.Ruler.MinorTick = 'on';
title('Lifted and Restricted Points Colored by Distance from Original Embedding', 'fontsize', 16);
xlabel('\psi_1', 'fontsize', 18);
ylabel('\psi_2', 'fontsize', 18);
ylabel('\psi_2', 'fontsize', 12);
