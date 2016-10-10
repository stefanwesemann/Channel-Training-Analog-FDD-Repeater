# This script is used to generate Figure 4, in the article:
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
library(reshape2)


# Load/define Simulator Functions ----------------------------------------------

if(!exists("simSsRSME", mode="function")) source("fun_defs.R")

simPwrIt <- function(M,rhoUdB,rhoDdB,R,TMrat,eps=0,maxNumIt=100){
    # Fix RNG seed for reproducable results
    set.seed(123)
    
    # Calculate pilot signal length
    T <- M*TMrat
    
    # Compute linear SINRs
    rhoU <- 10^(rhoUdB/10)
    rhoD <- 10^(rhoDdB/10)
    
    # Placeholder for storing simulation results
    numIt <- matrix(0, nrow=R)
    
    # Loop over all Monte Carlo realizations
    for (realIdx in c(1:R)){
        
        # Generate random (unitary) pilot matrix
        Phi <- matrix(1/sqrt(2)*(rnorm(T*M)+1i*rnorm(T*M)),ncol=M)
        Phi <- qr.Q(qr(Phi))
        
        # Generate UL & DL channels as well as noise realizations
        g <- matrix(1/sqrt(2)*(rnorm(M)+1i*rnorm(M)),ncol=1)
        
        h <- matrix(1/sqrt(2)*(rnorm(M)+1i*rnorm(M)),ncol=1)
        
        w <- matrix(1/sqrt(2)*(rnorm(T)+1i*rnorm(T)),ncol =1)
        
        N <- matrix(1/sqrt(2)*(rnorm(M*T)+1i*rnorm(M*T)),ncol =T)
        
        # Construct repeater Rx signal
        x <- sqrt(rhoD*T/M)*Phi%*%h+w
        
        # Construct FH Rx signal
        Y <- sqrt(rhoU/(rhoD+1))*g%*%Conj(t(x))+N
        
        # Pilot-reverse modulation step
        Z <- Y%*%Phi
        
        # Power iteration based estimator
        tmp<- rank1SVD(Z,N=maxNumIt,eps=eps)
        numIt[realIdx] <- tmp$numIt
    }
    
    return(numIt)
}


# Simulation parameters --------------------------------------------------------

R <- 1e3                     # Number of Monte Carlo realizations
TMrat <-1                    # M/tau ratio (must be an integer)
Mlist <- c(4,8,16,32,64,128) # List of FH antenna number
rhoUlist <- c(-10,0,10)      # List of simulated UL SINRs(dB), defined in Eq.(4)
epsList <- c(1e-1,1e-2)      # List of threshold parameters for power iterations
rhoDdB <- 100                # DL SINR (dB), defined in Eq.(2)


# Run Simulations --------------------------------------------------------------

# Placeholder for storing simulation results
data <- matrix(0,nrow=length(Mlist)*length(rhoUlist)*length(epsList),ncol=6)

# Loop over all parameter values
for (eIdx in seq(1,length(epsList))){
    for (rIdx in seq(1,length(rhoUlist))){
        for (mIdx in seq(1,length(Mlist))){
        
            # Call simulator for current parameter set
            results <- simPwrIt(Mlist[mIdx],rhoUlist[rIdx],rhoDdB,R,TMrat,
                                eps=epsList[eIdx])
            
            # Compute CDF and specific percentiles
            tmp <- calcCDF(results,numBins=1e3,plotFlag=0)
            
            rowIdx <- mIdx+(rIdx-1)*length(Mlist) + 
                           (eIdx-1)*length(Mlist)*length(rhoUlist)
            
            data[rowIdx,1] <- Mlist[mIdx]
            data[rowIdx,2] <- rhoUlist[rIdx]
            data[rowIdx,3] <- epsList[eIdx]
            data[rowIdx,4] <- tmp$xm
            data[rowIdx,5] <- tmp$x005
            data[rowIdx,6] <- tmp$x095
        
        }
    }
}

# create data frame for plotting
dfG <- data.frame(data)
colnames(dfG) <- c("x","r","e","y","lower","upper")

dfG$r <- factor(dfG$r,labels=c(
    expression(paste(rho[U]==-10,"dB")),
    expression(paste(rho[U]==0,"dB")),
    expression(paste(rho[U]==10,"dB"))))
dfG$e <- factor(dfG$e,labels=c(expression(delta==10^-2),expression(delta==10^-1)))


# Plot Simulation Results ------------------------------------------------------
plotFileFlag <- 0   # flag for file output (as *.pdf)

if (plotFileFlag){
    scaleFactor <- 0.9
    pdf("numIt.pdf", width=scaleFactor*8, height=scaleFactor*6)
}

p <- ggplot(dfG, aes(x=log2(x),y=y)) + theme_bw() + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha = 0.3,fill="#56B4E9")+
    geom_line(colour="#56B4E9",size=0.5)+
    scale_x_continuous(breaks = log2(Mlist),labels = Mlist)+
    scale_y_continuous()+
    ylab("Number of Iterations")+
    xlab("M")
    
p <- p + theme(legend.title=element_blank(),
               legend.key = element_rect(colour = NA),
               legend.text.align = 0,
               text = element_text(size=16),
               panel.grid.minor = element_blank())

p <- p + facet_grid(e ~ r, labeller = label_parsed) +
     theme(strip.text = element_text(size=16),
           strip.background = element_rect(fill="white", colour="white",size=1),
           axis.text.x = element_text(size=10),
           axis.text.y = element_text(size=10))

p

if (plotFileFlag){
    dev.off()
}