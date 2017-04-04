# This script is used to generate the left side of Figure 3, in the article:
#
# Stefan Wesemann, Thomas L. Marzetta, "Channel Training for Analog FDD 
# Repeaters: Optimal Estimators and Cramér-Rao Bounds"
# Download article: https://arxiv.org/pdf/1610.03260v2.pdf
#
# This is version 1.0 (Last edited: 2017-04-04)
#
# License: This code is licensed under the GPLv2 license. If you in any way
# use this code for research that results in publications, please cite our
# original article listed above.


# load required libraries ------------------------------------------------------
# Please make sure that the corresponding packages are installed!

rm(list=ls())

# Install function for packages    
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

packages(ThreeWay)
packages(latex2exp)
packages(pracma)
packages(plyr)
packages(ggplot2)


# Load Simulator Functions -----------------------------------------------------

if(!exists("simSsRSME", mode="function")) source("fun_defs.R")


# Simulation parameters --------------------------------------------------------

R <- 1e3                     # Number of Monte Carlo realizations
TMrat <-1                    # M/tau ratio (must be an integer)
rhoDLdB <- 20                # DL SINR (dB), defined in Eq.(2)
eps <- 1e-1                  # Threshold parameter for power iterations
MList <- c(4,16,64)          # List of FH antenna number
rhoUList <- seq(-10,30,by=5) # List of simulated UL SINRs(dB), defined in Eq.(4)
calcModeList <- c("ULCRB","ULSVD","ULPWR")
                             # List of simulated methods; i.e., 
                             # Eq.(25), Eq.(20) and Algorithm (1)


# Run Simulations --------------------------------------------------------------

# Create parameter grid (as a data frame for efficient plot generation)
par.grid <- expand.grid(MList,rhoUList,calcModeList)
colnames(par.grid)<- list("M","rhoULdB","calcMode")

# Call simulator for each grid point
results <- mdply(par.grid,simSsRSME,R=R,rhoDLdB=rhoDLdB,TMrat=TMrat,eps=eps)

# Re-cast parameters
results$M <- factor(results$M)
results$rhoU<- factor(results$rhoU)
results$calcMode<- factor(results$calcMode)
colnames(results)[colnames(results) == 'V1'] <- "RMSE" 


# Plot Simulation Results ------------------------------------------------------
plotFileFlag <- 0     # flag for file output (as *.pdf)

if (plotFileFlag){
    scaleFactor <- 0.9
    pdf("UL_rmse.pdf", width=scaleFactor*8, height=scaleFactor*6)
}

p <- ggplot(results,
            aes(x=rhoU,
                y=RMSE,
                colour=calcMode,
                shape=calcMode,
                linetype=M,
                group=interaction(M,calcMode)))+
    geom_line(size=0.5)+
    geom_point(size=4)

p <- p + theme_bw()

p <- p + scale_y_log10() +
    scale_x_discrete(breaks=seq(-10,30,by=10)) +
    annotation_logticks(base = 10,sides="lr")

p <- p + xlab(expression(paste(rho[U]," (dB)")))+ ylab("RMSE (rad)")

p <- p + theme(legend.title=element_blank(),
               legend.key = element_rect(colour = NA),
               legend.text.align = 0,
               text = element_text(size=16),
               panel.grid.minor = element_blank())

# Change the legend
p <- p + scale_linetype_discrete(name  =NULL,
                                 breaks=c(4,16,64),
                                 labels=c("M=4", "M=16","M=64")) +
    scale_shape_discrete(name  =NULL,
                         breaks=c("ULCRB","ULSVD","ULPWR"),
                         labels=c("Cramér Rao Bound (25)", 
                                  "Simulation: SVD (LAPACK)",
                                  expression(paste(
                                      "Simulation: Power Iteration, ",
                                      delta==10^-1))),
                         solid = FALSE) +
    scale_colour_discrete(name  =NULL,
                          breaks=c("ULCRB","ULSVD","ULPWR"),
                          labels=c("Cramér Rao Bound (25)", 
                                   "Simulation: SVD (LAPACK)",
                                   expression(paste(
                                       "Simulation: Power Iteration, ",
                                       delta==10^-1))))

p <- p + theme(legend.box.just = "left",
               legend.justification=c(0,0), 
               legend.position=c(0,0))

p <- p +guides(colour = guide_legend(order = 2),
               shape = guide_legend(order = 2),
               linetype = guide_legend(order = 1))

p

if (plotFileFlag){
    dev.off()
}
