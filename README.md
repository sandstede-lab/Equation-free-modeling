# eq-free
Enabling equation-free modeling via diffusion maps by Tracy Chin, Jacob Ruth, Clayton Sanford, Rebecca Santorella, Paul Carter, and Bjorn Sandstede

Generate traffic data: genTrafficData.m 

Compute a diffusion map:

Computing the bifurcation diagram
1. Using the microsystem: microBifurcation.m
    
2. Using the standard deviation as a marco variable: eqFreeBifurcation.m 

3. Using a 1D diffusion map applied to phase-shifted traffic profiles: eqFreeDiffBifurcation.m 

4. Using a 2D diffusion map applied to traffic profiles: eqFreeDiffBifurcation2D.m 


Data files
1. microBif.mat - Matlab datafile containing the data points on the bifurcation curve by sigma and v0 saved as 'bif'

2. refStates.mat - Matlab datafile containing the two starting references states 'ref1' and 'ref2' as well as their velocities 'v0_base1' and 'v0_base2'

