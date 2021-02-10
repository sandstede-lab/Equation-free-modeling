%% 1D diffusion map



%% 2D diffusion map


[~, max1] = max(allData,[],1);  % locate the max headway for each data point
% plot eigenvector 1 vs eigenvector 2, colored by max headway location
figure(1);
scatter(evecs(:,1), evecs(:,2), 100, max1,'.'); hold on;
color = colorbar;
xlabel(color, 'Wave Position', 'fontsize', 14)
colormap(jet);
xlabel('\psi_1', 'FontSize',14);
ylabel('\psi_2', 'FontSize',14);
title('\psi_1 vs. \psi_2 Colored by Locations of the Max Headways','FontSize',14);
drawnow;
% plot eigenvector 1 vs eigenvector 2, colored by standard deviation
figure;
scatter(evecs(:,1), evecs(:,2), 100,  std(allData),'.');
colorbar;
color = colorbar;
xlabel(color, '\sigma', 'fontsize', 14)
colormap(jet);
xlabel('\psi_1', 'FontSize',14);
ylabel('\psi_2', 'FontSize',14);
title('\psi_1 vs. \psi_2 Colored by Standard Deviation of the Headways','FontSize',14);
