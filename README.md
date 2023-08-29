# Stress_Patterning_and_Finite_Time_Elastic_Singularities_in_a_Renewable_Active_Elastomer
This contains codes for the manuscript 'Stress Patterning and Finite Time Elastic Singularities in a Renewable Active Elastomer'

The file non_linear_eqn.m is written in MATLAB and the files one_species.py, two_species.py is written in Python using the package Dedalus.

1) The file non_linear_eqn.m is used to solve the system of non-linear partial differential equations (1: as in the main text) and generate the numerical phase diagram. The file can be compiled and run using MATLAB and a movie in .avi format is generated upon completion of run-time. The code can be run in different parameter regimes (as shown in the numerical phase diagram) to generate the different phases. Movies for the travelling wave, swap phase and temporal coexistance of travelling wave and swap phase in a binary mixture of stresslets have been generated using this code i.e. Movies S4, S5 and S6.

2) The profiles for the segregation regime is analysed using codes written in Python using the pseudo-spectral solver Dedalus. Snapshots (containing the values of the variables and quantities of interest) are saved in HDF5 files. The file 'one_species.py' is for a single stresslet on an elastomer. This code can be run for the parameter regimes (either the segregation or the wave regime) mentioned in the Supplementary Material. The file 'two_species.py' is for two stresslets on an elastomer. Movie S1, S2, S3, S7 and S8 are generated using data generated from these codes.


3) 

