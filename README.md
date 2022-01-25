# Taka_A1_A2_only
Comparison A1 and A2 progenitor cells

This R project contains the code to generate Figures for the publication titled 

###Specification of CNS macrophage subsets occurs postnatally in defined anatomical niches

The project can be started in RStudio by executing the file Taka_A1_A2_only.Rproj. The folder structure is as follows:
- **bin/** contains the R scripts to run the project. The number in front of the file name indicates the order in which the code should be executed, i.e. start with 0_setup.R to setup the folder structure and install required packages, then 1_load_counts.R, 2_seurat.R, and 4_plotting.R. functions.R contains custom functions for plotting and colors used for the plots.

- **data/** contains the data for the project. The counts files are contained in the **counts** folder. All the data associated with the project is saved and loaded from here. The seurat file is called all.RData.

- **plots/** contains the plots used in the associated publication.

With regards to data exclusion, the 2 main points are noteworthy. The counts data contains cells from the so called EMP gate. These are erythromyeloid progenitors that were acquired and analyzed with the data, but not shown in the paper as the final focus was on A1 and A2 cells. Furthermore, some cells appeared to be damaged. These cells are documented in the file data/damaged-cells.csv and were excluded before running the Seurat algorithm.

The code is deposited on github under: https://github.com/rsankowski/Taka_A1_A2_only

In case of questions, please reach out to me at any time under roman.sankowski@uniklinik-freiburg.de