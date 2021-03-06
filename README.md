# Enabling Equation-Free Modeling via Diffusion Maps 

Equation-free modeling estimates macro-level behavior through a coarse time stepper in three steps: (1) lift: build the microstate from the macrostate, (2), evolve: simulate the microsystem for short bursts, and (3) restrict: estimate the macrostate from the evolved microstate. 

In this work, we use diffuson maps to identify macroscopic variables in a traffic model and define appropriate lifting and restriction operators. We then use these operators to compute and continue traffic jam solutions. 

**Folder navivgation:**  
**1.) src** contains the source code  
**2.) data** contains raw data files in .csv files
**3.) results** contains .csv files of results 
 
## Generate traffic data
Data used in the diffusion maps is created with genTrafficData.m, and saved in the /data directory. 

## Compute a diffusion map
Diffusoin maps for the traffic data are computed with and explored through createDiffMaps.m.

## Test the lifting and restriction operators
We test the accuracy of the lifting and restriction operators with the script testOperators.m.

## Computing the bifurcation diagram
**1.)** Using the microsystem: microBifurcation.m
    
**2.)** Using the standard deviation as a marco variable: eqFreeBifurcation.m 

**3.)** Using a 1D diffusion map applied to phase-shifted traffic profiles: eqFreeDiffBifurcation.m

**4.)** Using a 2D diffusion map applied to traffic profiles: eqFreeDiffBifurcation2D.m 

**BibTex Citation:**  
```
@article {eq-free-traffic,  
	author = {Chin, Tracy and Ruth, Jacob and Sanford, Clayton and Santorella, Rebecca and Carter, Paul and Sandstede, Bjorn},  
	title = {Enabling Equation-Free Modeling via Diffusion Maps},  
	year = {2021},  
}
```

