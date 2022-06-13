# Active_force_patterning
This contains codes for the manuscript 'Emergence of Tension Chains and Active Force Patterning'

The file non_linear_eqn.m is written in MATLAB and the file stresslets_segregation.py is written using Dedalus and Python. Sections of the code are commented and details on how to compile the code are mentioned in the files.

1) The file non_linear_eqn.m is used to solve the system of non-linear partial differential equations (1: as in the main text) and generate the numerical phase diagram. The file can be compiled and run using MATLAB and a movie in .avi format is generated upon completion of run-time. The code can be run in different parameter regimes (as shown in the numerical phase diagram) to generate the different phases. Movies for the travelling wave, swap phase and temporal coexistance of travelling wave and swap phase in a binary mixture of stresslets have been generated using this code i.e. Movies S3, S4, S5.

2) The segregation regime is analysed using codes written in the Dedalus scheme. Snapshots (images) in time are generated upon succesful completion of the run time and they can be stitched together to generate a movie using the ffmpeg package.
The code can be run for the parameter regimes mentioned in the Supplementary Material. Movies S1, S2, S6, S7, S8, S9 have been generated using this code. 
