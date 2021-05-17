# pyPROFILE
Python implementation of the basic functionality of the PROFILE software for biogechemistry by Berg et al., 1998

The code solves the diffusion equation for the concentration and returns the optimal profile of net rate of production R that best fits an input experimental profile. This corresponds to the procedure described by eq. (1) to (10) of Berg et al., 1998.
Unlike Berg et al., this code doesn't optimize for the number of zones. However, I don't think it is necessary features because (1) modern computers can compute the equation for >1000 in a matter of seconds; (2) dividing the profile into a few zones where R is piecewise constant can hide the fact that R varies continuously. If users resquest this function anyway, I'm open to implementing it later. 

## Required python libraries
All the required libraries are part of the standard scientific python stack. I recommennd installing [Anaconda](https://www.anaconda.com/). Required libraries:
- numpy
- scipy
- pandas
- matplotlib

## Example

The file `Example.ipynb` is a notebook that demonstrates how to call the optimization function and how to extract/visualize the results. The file `./Data/Example.csv` is an example of input file. Input files can contain several concentration profiles sampled at the same depth (i.e., replicas of experiment).


## Boundary conditions
For the moment, the only available condition is fixed concentration (i.e. Dirichlet) at top and bottom. The values are taken as the first and last value of the data. Note that flow boundary conditions (i.e. Neumann are implemented but not yet seriously tested).


References: 
Berg, P., Risgaard-Petersen, N., Rysgaard, S., 1998, Interpretation of measured concentration profile in sediment pore water, Limnol. Oceanogr., 43(7), 1500-1510.
