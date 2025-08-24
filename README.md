## Evaporation of Droplets - Numerical Model

### Description

This repository contains a MATLAB implementation of a numerical model that simulates the evaporation of multicomponent spherical droplets with variable physical properties (including density). The model solves the mass and energy balance equations iteratively, updating properties at each time step and in each part of the domain.

### Requirements

- MATLAB R2023a, including Optimization Toolbox

### Usage

1. Clone the repository
2. Open MATLAB and set the working directory to the project folder.
3. Adjust parameters in the main script: `Run_droplet_test.m`
4. Adjust parameters in: `test/Droplet/droplet_test.m`
5. Run the main script.
6. Adjust parameters in the script: `test/TransformSolution.m`
7. Run the script.

### Folder Structure

`src`                  - Core model functions \
`test`                 - Main simulation codes called from the initialization script and post-processing functions \
`results`              - Folder where the results are saved (***) \
`Run_droplet_test.m`  - Main script to run the simulation \

(***) It is recommended, once the simulation is finished, to store the files generated in the results folder in another folder that is not included in the path to avoid errors in future simulations.

### Saved Files

The simulation generates three main files:
- `Model.mat`: Contains characteristic input data of the simulation before it starts, such as the number of components, type of fuels, polynomial order, etc.
- `Saved_vars.mat`: Stores selected relevant variables for each time step during the simulation.
- `Saved_solutions.mat`: Stores the complete simulation results at selected time intervals, as determined by the parameters defined in the initialization script. These results require further processing for visualization.

At the end of a successful simulation, the last two files are automatically saved as:
- `Saved_varsEnd.mat`
- `Saved_solutionsEnd.mat`

However, if the simulation is interrupted (due to an error or excessive run time), it can be resumed in parts from the last saved time step. To enable this, automatic backup files are generated throughout the simulation (the frequency of these backups is set in the initialization settings).
When processing the results of a simulation that was run in parts, the corresponding backup files must be merged.

For example, if the simulation created: `Saved_solutions20.mat`, `Saved_solutions40.mat`, `Saved_solutions60.mat`, `Saved_solutions80.mat`, `Saved_solutions100.mat`, and `Saved_solutionsEnd.mat` 
and the simulation was restarted from `Saved_solutions40.mat`, then only the following files should be combined: `Saved_solutions40.mat` and `Saved_solutionsEnd.mat`.
Similarly for: `Saved_vars40.mat` and `Saved_varsEnd.mat`

### Authors

Daniel Betrán (Fluid and Energy Engineering Laboratory, LIFEn) \
Álvaro Gutiérrez (Fluid and Energy Engineering Laboratory, LIFEn) \
Manuel Colera (Universidad Politécnica de Madrid) \
Álvaro Muelas (Fluid and Energy Engineering Laboratory, LIFEn) \
Javier Ballester (Universidad de Zaragoza)

### Cite as

(Por determinar)


