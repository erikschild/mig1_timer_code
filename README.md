# mig1_timer_code
This repository contains all data and code necessary to reproduce mathematical modeling, statistical testing and figures for [[citation manuscript]]

Matlab code:

To plot Fig 2D, 3A:
- Run fig23_process1.m: fits activator model to control data (fig234_data_ctrl.txt)
- Run fig23_process2.m: fits repressor model to control data (fig234_data_ctrl.txt)
- Run fig23_plot: plots figure panels

To plot Fig 4E,F:
- Run fig4_process1.m: finds division time from control data (fig234_data_ctrl.txt)
- Run fig4_process2.m: finds mean size ratio from wild-type data (fig4_data_ratio.txt)
- Run fig4_process3.m: fits activator model with division to control data (fig234_data_ctrl.txt)
- Run fig4_plot: plots figure panels

To compute statistics for Fig 5B:
- Run fig5_process: runs ANOVA on data (fig5_data_ctrl.txt, fig5_data_gof.txt, fig5_data_lof.txt)

R code:

To plot figure 1B+C+E, 2B+C, 3C, 4A+B+C+D, 5A+B+S1C, S1A+B, and calculate corresponding statistics:
- Run corresponding .R file
- for Figure 2B+C.R and Figure 3C.R, note that the greek uppercase delta symbol may not properly show in the code, depending on encoding. UTF-8 encoding allows for greek letters
