# Deep scoping: a future-proof breeding strategy
## Running the code
- Install the package hypred (zip file included)
- Install the package GSSimTPUpdate (zip file included)
- Open the R project OracleSelection.Rproj
- Run Create_directory.R
- Run MakeGenome_File.R
- Run Make_Base_Population.R
- Run Make_Populationlist.R
- Run one of the 'run_experiment' files in the main directory
- Transfer the simulation results from the 'own_results' directory to the 'data' directory
- Run the correct 'make_figure' script

## data 
This directory contains the original base population used for each simulation study and several directories where the simulation results can be stored.

## make_figures
This directory contains the R script used to visualize the results. The simulated data should be put in the data/\<method\> directory. 

## own_results
Empty directory where simulation results will be saved. This directory is creating by running the script Create_directory.R.

## MakeGenome_File 
Makes a list of n different genomes that can be used in the run_experiment files. At each iteration, each method will use the same genome, making it possible to compare the different methods.

## Genome
Directory containing an example files of genomes that can be used in the run_experiment files.

## Make_Base_Population
Simulates truncation selection for 5, 10, 15, and 20 breeding cycles and saves the different populations into the Population directory.

## Make_Populationlist
Imports the population that was created with Make_Base_Population.R and save the different simulation parameters in a list that can be used in the run_experiment files.

## Population
A directory containing the base populations after 5, 10, 15 and 20 breeding cycles of truncation selection. Each run_experiment file requires a population. To create a population first run the Make_Base_Population.R followed by the Make_Populationlist.R.

## run_experiment files
R scripts used to simulate a population following the deep scoping method, HUC method with bridging, scoping method, and population merit method. The results will be saved in the directory 'own_results'. 
The user will have to load the correct genome from the genome directory and choose the number of nodes that are available for calculation. 

## R
Directory containing functions that are used in the run_experiment files.

## Run files
#run_experiment_Scoping_Method.R
Main file to simulate a breeding population using the scoping method. This file can be used to generate the different TP update methods as well as the scoping method that was compared with the oracle parental selection.

#run_experiment_Truncation_Selection.R
Main file to simulate a breeding population using truncation selection. This file can be used to generate the different TP update methods as well as the scoping method that was compared with the oracle parental selection.

#run_experiment__OGM
Run the OGM method 

#run_experiment_Chimeric_Scoping
Run Chimeric Scoping

#run_experiment_DeepScoping_method
Run Deep Scoping

#run_experiment_Oracle_Parental_Selection
Run Oracle parental selection

#run_experiment_Oracle_TP_Update_Scoping_Method
Run oracle TP update method using the scoping method

#run_experiment_Oracle_TP_Update_Truncation_Selection
Run oracle TP update method using truncation selection