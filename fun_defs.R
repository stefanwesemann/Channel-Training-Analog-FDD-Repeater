# This script contains function definitions that are needed for generating
# Figure 2, in the article:
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


# Define binary operator for calculating the subspace distance
'%A%' <- function(a,b) {
    a<- matrix(a,ncol=1)
    b<- matrix(b,ncol=1)
    return(acos(pmin(1,Mod(Conj(t(a))%*%b)/(norm(a,"2")*norm(b,"2")))))
}


# Function calculates the inverse Fisher information metric
invFIM <- function(M,T,gamma){
    return((M-1)*(1+gamma)/(T*gamma^2))
}


# Algorithm (1): Rank-1 SVD based on power iterations
rank1SVD <- function(Z,eps=1e-2,N=1e2){
    # Find column with largest norm
    strongestSmpIdx <- which.max(colMeans(Mod(Z)^2))
    
    # Initialize g
    gbar <- Z[,strongestSmpIdx]
    gEst <- gbar/norm(gbar,"2")
    
    # Iterative updates for the estimates of g and h
    itVec <- matrix(0,nrow=10,ncol=1)
    pwrItCntr <- 0
    for (itIdx in seq(1,N)){
        pwrItCntr <- pwrItCntr +1
        
        hbar <- Conj(t(Z))%*%gEst
        hEst <- hbar/norm(hbar,"2")
        
        gbar <- Z%*%hEst
        gEstOld <- gEst
        gEst <- gbar/norm(gbar,"2")
        
        # Check stopping criterion
        if ((gEstOld %A% gEst) <=eps) break
    }
    return(list(gEst=gEst,hEst=hEst,numIt = pwrItCntr))
}


# Perform Monte Carlo simulations 
simUlDl <- function(M,T,P,R,calcMode,beta,alpha,eps=1e-2){
    
    # Fix RNG seed for reproducable results
    set.seed(123)
    
    # Placeholder for storing simulation results
    dssUL <- matrix(0,nrow=R, ncol=1)
    dssDL <- matrix(0,nrow=R, ncol=1)

    # Loop over all Monte Carlo realizations
    for (realIdx in c(1:R)){
        # Generate random (unitary) pilot matrix
        Phi <- matrix(1/sqrt(2)*(rnorm(T*M)+1i*rnorm(T*M)),ncol=M)
        Phi <- qr.Q(qr(Phi))
        
        # Generate UL & DL channels as well as noise realizations
        g <- matrix(sqrt(beta)/sqrt(2)*(rnorm(M)+1i*rnorm(M)),ncol=1)
        
        h <- matrix(sqrt(beta)/sqrt(2)*(rnorm(M)+1i*rnorm(M)),ncol=1)
        
        w <- matrix(1/sqrt(2)*(rnorm(T)+1i*rnorm(T)),ncol =1)
        
        N <- matrix(1/sqrt(2)*(rnorm(M*T)+1i*rnorm(M*T)),ncol =T)
        
        # Construct repeater Rx signal
        x <- sqrt(T/M*P)*Phi%*%h+w
        
        # Construct FH Rx signal
        Y <- sqrt(alpha)*g%*%Conj(t(x))+N
        
        # Pilot-reverse modulation step
        Ytilde <- Y%*%Phi
        
        # ML subspace estimation
        if ((calcMode=="ULSVD")||(calcMode=="DLSVD")){
            # SVD-based estimator
            USV <- svd(Ytilde,nu=1,nv=1)
            gEst <- USV$u[,1]
            hEst <- USV$v[,1]
        } else if ((calcMode=="ULPWR")||(calcMode=="DLPWR")){
            # Power iteration based estimator
            tmp <- rank1SVD(Ytilde,eps)
            gEst <- tmp$gEst
            hEst <- tmp$hEst
        }
        
        # Store simulation results
        dssUL[realIdx] <- gEst%A%g
        dssDL[realIdx] <- hEst%A%h
    }
    
    # Calculate root mean square error
    if ((calcMode=="ULSVD")||(calcMode=="ULPWR")){
        rmse <- sqrt(mean(dssUL^2))
    } else if ((calcMode=="DLSVD")|| (calcMode=="DLPWR")){
        rmse <- sqrt(mean(dssDL^2))
    } 
    return(rmse)
}


# Wrapper function for switching between Monte Carlo sim. and CRB computation
simSsRSME <- function(M,rhoULdB,calcMode,R,rhoDLdB,TMrat,eps=1e-2){
    
    # Calculate pilot signal length
    T <- M*TMrat
    
    # Compute linear SINRs
    rhoUL <- 10^(rhoULdB/10)
    rhoDL <- 10^(rhoDLdB/10)
    
    # Assume unit transmit power for FH (without loss of generality)
    P <- 1
    
    # Derive variance of the channel coefficients 
    beta <- rhoDL/P
    
    # Compute scaling factor for the repeater feedback signal 
    alpha <- rhoUL/(beta*(beta*P+1))
    
    # Compute effective SINRs
    rhoDLTilde <- T/M*beta*P
    rhoULTilde <- alpha*beta*(T/M*beta*P+1)
    
    
    if ((calcMode=="ULSVD")||(calcMode=="DLSVD")||
        (calcMode=="ULPWR")||(calcMode=="DLPWR")){
        # Perform Monte Carlo simulation
        rmse <- simUlDl(M,T,P,R,calcMode,beta,alpha,eps)
        
    } else if (calcMode=="ULCRB"){
        # Compute UL CRB based on Eq.(25)
        rmse <- sqrt(invFIM(M,M,M*rhoULTilde))
        
    } else if (calcMode=="DLCRB"){
        # Compute DL CRB based on Eq.(41)
        rmse <- sqrt(invFIM(M,1,M*rhoDLTilde)+
                         invFIM(M,M,M*rhoULTilde))
    } 
    
    return(rmse)
    
}



# Function calculates the log-likelihood of the channel gain (likelihood given
# by Eq.(56))
calcLLH <- function(ccIn,M,lambda,Qt){
    cc <- ccIn*Qt
    ss <- 0
    aVec <- matrix(0,nrow=M,ncol=1)
    for (ii in seq(1,M)){
        pp <- 1
        for (jj in seq(1,M)){
            if (ii!=jj){
                pp <- pp * cc/(cc+1)*(lambda[ii]-lambda[jj])
            }
        }
        aVec[ii] <- (exp(cc/(cc+1)*lambda[ii])-exp(cc/(cc+1)*lambda[1]))/pp
        ss <- ss + (exp(cc/(cc+1)*lambda[ii])-exp(cc/(cc+1)*lambda[1]))/pp
    }
    likelihood <- -M*log(cc+1) + log(max(c(1e-16,ss)))
    return(likelihood)
}



# Golden Section search method as given by Appendix C.3 in D. Bertsekas, 
# "Nonlinear Programming", Athena Scientific, 1995.
goldSecSearch <- function(M,lambda,Qt,lowerBound,upperBound,precision){
    
    # Define the function that needs to be minimized
    f <- function(x){return(-calcLLH(x,M,lambda,Qt))}
    
    # Define Golden Ratio
    tau <- (3-sqrt(5))/2
    
    # Initialize bounds
    a  <- lowerBound
    ab <- upperBound
    
    while(1) {
        
        b  <- a + tau*(ab-a)
        bb <- ab - tau*(ab-a)
        
        ga <- f(a)
        gab <- f(ab)
        
        if (is.nan(gab)) {
            return(NaN)
            break
        }
        
        gb <- f(b)
        gbb <- f(bb)
        
        if (gb < gbb){
            if (ga <= gb){
                a <- a
                ab <- b
            } else {
                a <- a
                ab <- bb
            }
        } else if (gb>gbb){
            if (gbb>gab) {
                a <- bb
                ab <- ab
            } else {
                a <- b
                ab <- ab
            }
            
        } else {
            a <- b
            ab <- bb
        }
        
        if (Mod(a-ab)< precision) {break}
        
    }
    return(1/2*(a+ab))
}



# Simplified Monte Carlo simulation chain for UL channel gain estimation
simULChanGain <- function(calcMode,rhoUdB,M,R){
    
    # Fix RNG seed for reproducable results
    set.seed(123)
    
    # Compute linear UL SINR
    rhoU <- 10^(rhoUdB/10)
    
    # Assume unit (effective) transmit power for repeater
    Qt <- 1
    
    # Derive variance of the channel coefficients 
    beta <- rhoU/Qt
    
    # Placeholder for storing simulation results
    ulChanGain <- matrix(0,nrow=R, ncol=1)
    
    # Loop over all Monte Carlo realizations
    for (realIdx in c(1:R)){
        # Construct effective repeater Tx signal
        xtilde <- matrix(sqrt(Qt)/sqrt(2)*(rnorm(M)+1i*rnorm(M)),ncol=1)
        
        # Generate UL channel vector and Perturbation realizations
        g <- matrix(sqrt(beta)/sqrt(2)*(rnorm(M)+1i*rnorm(M)),ncol=1)
        N <- matrix(1/sqrt(2)*(rnorm(M^2)+1i*rnorm(M^2)),ncol =M)
        
        # Construct effective FH receive signal
        Ytilde <- g%*%Conj(t(xtilde))+N
        
        # Compute outer-product of received signal and determine the eigenvalues
        YtYtH <- Ytilde%*%Conj(t(Ytilde))
        YtYtHevd <- eigen(YtYtH)
        lambda <- YtYtHevd$values
        
        # Switch between ML estimator and SCM-based estimator
        if (length(grep("ML",calcMode))>0){
            chi <-goldSecSearch(M,lambda,Qt,1e-2,3e2,1e-4)
        } else if (length(grep("SCM",calcMode))>0){
            chi <- (lambda[1]/M-1)/Qt
        }
        
        # Store the result
        ulChanGain[realIdx] <- chi/M
    }
    
    # Switch between the computation of the first or second order statictic 
    if (length(grep("MEAN",calcMode))>0){
        # Check for numerial problems
        if (sum(is.nan(ulChanGain))/R<=0){
            # Calculate relative bias
            return((mean(ulChanGain, na.rm = TRUE)-beta)/beta)
        } else {
            return(NaN)
        }
        
    } else if (length(grep("VAR",calcMode))>0){
        # Check for numerial problems
        if (sum(is.nan(ulChanGain))/R<=0){
            # Calculate relative error variance
            return(1/beta^2*var(ulChanGain-mean(ulChanGain, na.rm = TRUE), na.rm = TRUE))
        } else {
            return(NaN)
        }
    }
}


# Calculate cumulative distribution function and 5th/50th/95th-percentile
calcCDF <-function(samples,numBins=1e2,plotFlag=0){
    
    samples <- matrix(samples,ncol=1)
    minVal <- 0.9*min(samples)
    maxVal <- 1.1*max(samples)
    
    bins <- seq(minVal,maxVal, by=((maxVal-minVal)/numBins))
    w<-hist(samples,bins, plot=FALSE)
    xTicks <- w$mids
    counts <-w$counts
    pdf <- counts/sum(counts)
    cdf <- cumsum(pdf)
    
    if (plotFlag==TRUE) {
        plot(xTicks,cdf,type="l",col="blue",
             xlim=c(min(xTicks),max(xTicks)),
             ylim=c(0,1),
             xlab="X",
             ylab="CDF(X)")
        grid()
    }
    
    idx05 <- min(which(cdf>=0.5))
    x05 <- mean(xTicks[c(idx05-1,idx05)])
    
    idx005 <- min(which(cdf>=0.05))
    x005 <- mean(xTicks[c(idx005-1,idx005)])
    
    idx095 <- min(which(cdf>=0.95))
    x095 <- mean(xTicks[c(idx095-1,idx095)])
    
    xm <- mean(samples)
    
    return(list(xTicks=xTicks,pdf=pdf,cdf=cdf,x05=x05,x005=x005,x095=x095,xm=xm))
}