# This script is used to generate the left side of Figure 3, in the article:
#
# Stefan Wesemann, Thomas L. Marzetta, "Channel Training for Analog FDD 
# Repeaters: Optimal Estimators and Cramér-Rao Bounds"
# Download article: tbd
#
# This is version 1.0 (Last edited: 2016-10-10)
#
# License: This code is licensed under the GPLv2 license. If you in any way
# use this code for research that results in publications, please cite our
# original article listed above.

# load required libraries ------------------------------------------------------
# Please make sure that the corresponding packages are installed!

rm(list=ls())

library(ThreeWay)
library(latex2exp)
library(pracma)
library(plyr)
library(ggplot2)


# Load Simulator Functions -----------------------------------------------------

if(!exists("simULChanGain", mode="function")) source("fun_defs.R")


# Simulation parameters --------------------------------------------------------

R <- 1e3                     # Number of Monte Carlo realizations
rhoUList <- seq(-10,10,by=2) # List of simulated UL SINRs(dB), defined in Eq.(4)
MList <- c(4,8)              # List of FH antenna number
calcModeList = c("ML_MEAN","SCM_MEAN")
                            # List of simulated methods; i.e., 
                            # ML estimator Eq.(52), SCM-based estimator Eq.(58) 


# Run Simulations --------------------------------------------------------------

# Create parameter grid (as a data frame for efficient plot generation)
par.grid <- expand.grid(MList,rhoUList,calcModeList)
colnames(par.grid)<- list("M","rhoU","calcMode")

# Call simulator for each grid point
results <- mdply(par.grid,simULChanGain,R=R)

# Re-cast parameters
results$M <- factor(results$M)
results$rhoU<- factor(results$rhoU)
results$calcMode<- factor(results$calcMode)
colnames(results)[colnames(results) == 'V1'] <- "beta" 


# Plot Simulation Results ------------------------------------------------------
plotFileFlag <- 0   # flag for file output (as *.pdf)

if (plotFileFlag){
    scaleFactor <- 0.9
    pdf("beta_mean.pdf", width=scaleFactor*8, height=scaleFactor*6)
}

p <- ggplot(results,
            aes(x=rhoU,
                y=(beta),
                colour=M,
                shape=calcMode,
                linetype=calcMode,
                group=interaction(calcMode,M)))+
    geom_line(size=0.5)+
    geom_point(size=4)

p <- p + theme_bw()

p <- p + xlab(expression(paste(beta," (dB)")))+ 
    ylab(expression(paste(beta^{-1},E(hat(beta))-1)))

p <- p + theme(legend.title=element_blank(),
               legend.key = element_rect(colour = NA),
               legend.text.align = 0,
               text = element_text(size=16),
               panel.grid.minor = element_blank())

# Change the legend
p <- p + scale_colour_discrete(name  =NULL,
                                 breaks=c(4,8),
                                 labels=c("M=4", "M=8")) +
    scale_linetype_discrete(name  =NULL,
                         breaks=c("ML_MEAN","SCM_MEAN"),
                         labels=c("ML Estimator (52)", 
                                  "SCM-based Estimator (58)")) +
    scale_shape_discrete(name  =NULL,
                         breaks=c("ML_MEAN","SCM_MEAN"),
                         labels=c("ML Estimator (52)", 
                                  "SCM-based Estimator (58)"),
                         solid = FALSE) 

p <- p + theme(legend.box.just = "right",
               legend.justification=c(0,0), 
               legend.position=c(0.5,0.65))

p <- p +guides(colour = guide_legend(order = 2),
               linetype = guide_legend(order = 1),
               shape = guide_legend(order = 1))

p

if (plotFileFlag){
    dev.off()
}