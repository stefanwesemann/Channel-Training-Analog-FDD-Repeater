Channel Training for Analog FDD Repeaters: Optimal Estimators and Cramér-Rao Bounds
==================

This is a code package is related to the following scientific paper:

Stefan Wesemann and Thomas L. Marzetta, “[Channel Training for Analog FDD Repeaters: Optimal Estimators and Cramér-Rao Bounds](https://arxiv.org/pdf/1610.03260v2.pdf),” submitted to IEEE Transactions on Signal Processing, October 2016, revised in April 2017

The package contains a simulation environment, based on R (https://www.r-project.org/), that reproduces all the numerical results and figures in the paper. 

## Abstract of Article

For frequency division duplex channels, a simple pilot loop-back procedure has been proposed that allows the estimation of the UL & DL channels at an antenna array without relying on any digital signal processing at the terminal side. For this scheme, we derive the maximum likelihood (ML) estimators for the UL & DL channel subspaces, formulate the corresponding Cramér-Rao bounds and show the asymptotic efficiency of both (SVD-based) estimators by means of Monte Carlo simulations. In addition, we illustrate how to compute the underlying (rank-1) SVD with quadratic time complexity by employing the power iteration method. To enable power control for the data transmission, knowledge of the channel gains is needed. Assuming that the UL & DL channels have on average the same gain, we formulate the ML estimator for the channel norm, and illustrate its robustness against strong noise by means of simulations.


## Content of Code Package

The article contains 5 simulation figures, numbered 3 (left & right), 4, and 5 (left & right). These are respectively generated by the R scripts figure_3_left.R, figure_3_right.R, figure_4.R, and figure_5_left.R, figure_5_right.R. The repository also contains the R script fun_defs.R which contains the definition of all common functions.

All scripts require the R packages "Threeway", "latex2exp", "pracma", "plyr", "reshape2" and "ggplot2". The code has been tested with RStudio and R version 3.2.1.

See each file for further documentation.


## Usage

Download the repository to a local directory, and set R's working directory to the same folder. It is important that all files are located in the same folder, otherwise the scripts will not find the underlying function definitions.


## Acknowledgements

This work has been conducted within Bell Labs’ FutureCell (FCell) project. The authors wish to thank project initiators T. Klein, A. Pascht and O. Blume. They are also grateful for the helpful comments of V. Suryaprakash.

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
