# o_lmg
Original Light Mediated Growth Simulation (No Device Physics). Development since December 2013.
Author Name: Nicolas A. Batara
Author Email: nick.batara@gmail.com

This set of scripts is designed to simulate photoelectrodeposition of thin semiconductor films. Simulaitons evolve through iterative process of calculating light absorption in FDTD and then updating structure based on monte carlo distribution function.

Supported Operating Systems:
Mac OSX
Linux

Software requirements: 
Lumerical FDTD
Matlab

Essential Files:

lmgSetup.m : Prompts user for simulation parameters and assembles simulation directories

lmg.fsp : Lumerical simulation file which includes modified analysis groups and optical constants for growth material

setup.lsf : Sets up initial simulation

nCores.txt : Specifies number of simulation cores to use

lmgAnalysis.m : Analyzes simulation results

run_lmg.m : Runs simulations in each directory by calling lmg.m

lmg.m : Main script which keeps track of simulation status, calls lmg_update.lsf and lmg_extract.lsf to 

lmg_update.lsf : Creates next .fsp iteration using the previous iteration and a structure file 

lmg_extract.lsf : Extracts absorption data from .fsp simulation and returns simulation to layout mode in order to save disk space.

