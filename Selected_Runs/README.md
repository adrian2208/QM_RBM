# Selected runs

This folder contains selected results from the program and how to reproduce them.								

## Figures

Shows Radial distribution, spatial distribution and energy as a function of ADAM learning step for a system of two interacting particles in two dimensions using importance sampling

## Procedure
1. In QM_RBM.cpp : 
   - Insert SelectedRuns(); into the main{} scope
   - Comment out any other functions in the main{} scope
   - Delete any potential files in the Output directory
   - compile from CMakeLists.txt in x64- Release mode and run QM_RBM.exe
2. In DataHandler.py :
   - At the mottom of the file, comment out any existing function calls
   - insert plot_Selected_Runs()
   - run the python program
   - the resulting graph should be saved to the Results folder, but not in the Report_Results folder
