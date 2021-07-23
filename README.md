# OLM_IS3control
"Deciphering how interneuron specific 3 cells control oriens lacunosum-moleculare cells to contribute to circuit function"
==============================================================================

This repository contains several self-contained folders for running the simulations described in "Deciphering how specialized interneuron-specific cell types contribute to circuit function". Accordingly, folders are labeled corresponding to the figure panels in the manuscript. Code was tested most recently in Python 3.7.6 on a Macbook Pro, MacOS Big Sur 11.2.2, and so some modification of these instructions may be needed for other versions and operating platforms.

-------------------------------------
Folders "Fig2A", "Fig2B", and "Fig2C"
-------------------------------------
Compile .mod files:
nrnivmodl

Choosing the OLM cell model:
In init_model.hoc adjust the "cell" variable (line 2) to numeric values of either 1 or 2, corresponding to cells 1 and 2 in the manuscript. Note that in the code, cell 1 is also referred to as "Gormadoc" and cell 2 is also referred to in "Isembard".

Run Simulations:
python init.py

Plotting:
Plots can be viewed mainly in the "PLOTfiles" folder, though FFT plots are saved separately in the "PLOTFFTfiles" folder.


-------------------
Folder "Fig3_and_4"
-------------------
Compile .mod files:
nrnivmodl

Choosing the OLM cell model:
In init_model.hoc adjust the "cell" variable (line 2) to numeric values of either 1 or 2, corresponding to cells 1 and 2 in the manuscript. Note that in the code, cell 1 is also referred to as "Gormadoc" and cell 2 is also referred to in "Isembard".

Run PRC simulations:
python init.py

Run code for plotting traces of interest:
python init_tracesofinterest.py

Run code for plotting cell 1 vs cell 2 (note that init.py needs to have been run for both cell 1 and 2 for this code to work):
python AnalyzeResulte.py

Plotting:
Plots for Fig 3 can be viewed mainly in the "PLOTfiles" folder, and plots for Fig 4 can be viewed in the "PLOTfiles_analysis" folder.

----------------------------------------------------------
Folders "Fig5ABCD", "Fig5EFGH", "Fig6ABCD", and "Fig6EFGH"
----------------------------------------------------------
Note that to run this code requires a high-performance computing (HPC) resources (e.g. we ran the code on NSG: https://www.nsgportal.org/). In light of this we have kept output results in these folders (i.e. saved in the "NPYFiles" folders) if a user would simply like to run analysis on the simulation output results and skip the simulation step.

Compile .mod files:
nrnivmodl

Choosing the OLM cell model:
No changes to the code is necessary here since the cell results are divided into separate master folders.

Run simulations (code is parallelized so requires HPC):
python init.py

Analyze and plot results:
python PlotResults.py

Plotting:
Plots can be found in the "Plots_" folders.

----------------------------------
Folders "Fig7ABCD", and "Fig7EFGH"
----------------------------------
Note that simulations for this code can be run on a local machine, but may take a while. In the simulation code there is commented-out code for running it in parallel however. We have kept output results from these simulations (i.e. saved in the "NPYFiles" folders) if a user would simply like to run analysis on the simulation output results and skip the simulation step.

Compile .mod files:
nrnivmodl

Choosing the OLM cell model:
For simulations, in init_model.hoc adjust the "cell" variable to numeric values of either 1 or 2, corresponding to cells 1 and 2 in the manuscript.
For analysis & plotting, in PlotResults.py, adjust the "Cell_Name" variable to either "Gormadoc" (cell 1) or "Isembard" (cell 2).

Run simulations (code is parallelized so requires HPC):
python init.py

Analyze and plot results:
python PlotResults.py

Plotting:
Plots can be found in the "Plots_" folders.
# OLM_IS3control
