The matlab scripts were written by Andrea Burke to allow calculation of the fraction of stratospheric sulfate in an ice core sulfate peak following the methods described in Burke et al., (2023). “High sensitivity of summer temperatures to stratospheric sulfur loading from volcanoes in the Northern Hemisphere”. PNAS. 

To calculate the fraction of stratospheric sulfate run the fstrat_solver.m file.

The additional files fstrat_MC.m and fstrat_equations.m must be in the same folder, as well as an excel data table with the sulfur isotope data formatted as in the Burke_2023_PNAS.xlsx data table. 

The name of the data table and the eruptions and cores wanted for the calculation need to be specified in the script before running.

Solutions are saved as a .mat file in the variable D.

These matlab files have been written and tested with Matlab version R2023a.
