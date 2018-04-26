# carbmud

These Matlab codes accompany the manuscript "The origin of carbonate mud" by Trower, Lamb, and Fischer.
The files were designed on Matlab R2017b on a Windows computer and do not need input to run; they can be called directly from the command line.

mud_fixedustar.m runs the mud production model at a fixed shear velocity (u_star) and plots the model with experimental data

mud_fixedD.m runs the mud production model at a fixed grain size (D) and plots the model with experimental data

mud.m runs the mud production model at varying grain size and shear velocity and plots the model with experimental data; this code may take a minute to run

susp_abrasion_calculations_mud.m is a set of calculations that accompanies the three preceding codes

settlingvelocity_mudvssand.m compares settling velocities and advection length scales for carbonate mud and sand

exptdata.mat contains the experimental data, which is used to plot against the model by the first three codes
