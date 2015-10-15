#-------------------------------------------------------------------------------
# Code to generate the numerous types of Q-Q plots 
#
# Adam Loy
# July 30, 2015
#-------------------------------------------------------------------------------

#######################
# DH Confidence bands #
#######################

sim_env <- function(x, conf = .95, line=FALSE){
  require(plyr)
  n <- length(x)
  P <- ppoints(x)[rank(x)]
  z <- qnorm(P)
  if (line) {
    require(HLMdiag)
    a <- as.numeric(HLMdiag:::qqlineInfo(x)[1])
    b <- as.numeric(HLMdiag:::qqlineInfo(x)[2])
  } else {
    a <- 0  # we know that the line should be the identity
    b <- 1
  }
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z)) * sqrt(P * (1 - P)/n)
  fit.value <- (a + b * z)
  upper <- (fit.value + zz * SE)
  lower <- (fit.value - zz * SE)
  RES <- data.frame(x = z, y = x, fit = fit.value, lower, upper)
  return(arrange(RES, x))
}

###############################################################
# The confidence band function                                #
# The function returns a list with the following:             #
# 1. lower,upper      the normal scale confidence bands       #
# 2. individual.alpha the final individual significance level #
###############################################################

QQ.cb<-function(x.sample,mu=0,sigma=1,M.sim=1000,alpha=0.05,plot=TRUE){
  
  n<-length(x.sample)
  upper.ci<-rep(NA,n)
  lower.ci<-rep(NA,n)
  p.value <-matrix(NA,nrow=n,ncol=M.sim)
  sim=NULL
  
  # simulate data
  for(i in 1:M.sim)
    sim<-cbind(sim,sort(runif(n)))
  
  # widen the CI to get a simultanoues 1-alpha CI
  for(i in 1:n){
    tmp<-pbeta(sim[i,],shape1=i,shape2=n+1-i)
    p.value[i,]<-apply(cbind(tmp,1-tmp),1,min)
  }
  
  critical.values<-apply(p.value,2,min)
  C.crit<-quantile(critical.values,prob=alpha)
  
  upper.ci<-qbeta(1-C.crit,shape1=1:n,shape2=n+1-(1:n))
  lower.ci<-qbeta(C.crit,shape1=1:n,shape2=n+1-(1:n))
  
  
  # now translate back to normal
  norm.upper<-qnorm(upper.ci)
  norm.lower<-qnorm(lower.ci)
  
  q.prob<-qnorm(ppoints(x.sample))
  z.sample<-(x.sample-mu)/sigma
  if(plot==TRUE){
    plot(q.prob,norm.upper,type="l",col="red",ylim=c(-3,3),xlim=c(-3,3),ylab="Sample Quantile",xlab="Sample Quantile")
    lines(q.prob,norm.lower,col="red")
    points(q.prob,sort(z.sample),pch=19,cex=0.6)
  }
  
  return(data.frame(x=q.prob, y=sort(z.sample), lower=norm.lower, upper=norm.upper))
}

###############################################################
# Q-Q plot functions                                          #
###############################################################

# Standard Q-Q plot
std_lineup <- function(dframe) {
  require(ggplot2)
  p <- ggplot(aes(x = x, y = y), data = dframe) + 
    geom_smooth(aes(x, fit), colour="grey50", se=FALSE, method="loess") +
    #  geom_abline(colour="grey50") +
    facet_wrap(~.sample, ncol=5) +
    geom_point() + 
    geom_ribbon(aes(x = x, ymin = lower.dh, ymax = upper.dh),alpha = .2)+
    labs(x = "", y = "") + theme_bw() + 
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"))
  p  
}


# Standard Q-Q plot with TS bands
std_ts_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  ggplot(aes(x =  x, y = y), data = dframe) + 
    #  geom_smooth(aes(naive1.qq.x, naive1.env.fit.value), colour="grey50", se=FALSE, method="loess") +
    geom_abline(colour = "grey50") +
    facet_wrap(~ .sample, ncol = 5) +
    geom_point() + 
    theme_bw() + 
    xlab("") + ylab("") +
    geom_ribbon(aes(x = x, ymin = lower.ts, ymax = upper.ts),alpha = .2)+
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"))
}


# Adjusted Detrended, DH bands
rot2_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  ggplot(aes(x = x, y = y - fit), data = dframe) + 
    geom_hline(yintercept = 0, colour = "grey30")+
    geom_point() + 
    geom_ribbon(aes(x = x, 
                    ymin = lower.dh - fit, 
                    ymax = upper.dh - fit), alpha = .2, data = dframe)+
    theme_bw() + 
    facet_wrap(~.sample, ncol=5) + 
    xlab("") + ylab("") +
    ylim(range(dframe$y)) +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"))
}


# Adjusted Detrended, TS bands
rot2_ts_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  ggplot(aes(x = x, y = y - fit), data = dframe) + 
    facet_wrap(~ .sample, ncol=5) + 
    geom_hline(yintercept=0, colour="grey30")+
    geom_point() + 
    geom_ribbon(aes(x = x, ymin = lower.ts - fit, ymax = upper.ts - fit), alpha = .2) +
    theme_bw() + 
    xlab("") + ylab("") +
    ylim(range(dframe$y)) +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"))
}


# Ordinary Detrended, DH bands 
rot_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  ggplot(aes(x = x, y = y - fit), data = dframe) + 
    facet_wrap(~ .sample, ncol = 5) + 
    geom_hline(yintercept = 0, colour = "grey30")+
    geom_point() + 
    geom_ribbon(aes(x = x, ymin = lower.dh - fit, ymax = upper.dh - fit), alpha = .2) +
    theme_bw() + 
    xlab("") + ylab("") +
    #            xlab("Normal Quantiles") + 
    #            ylab("Sample Quantiles") +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"))
}

# Ordinary Detrended, TS bands
rot_ts_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  ggplot(aes(x = x, y = y - fit), data = dframe) + 
    facet_wrap(~ .sample, ncol=5) + 
    geom_hline(yintercept=0, colour="grey30")+
    geom_point() + 
    geom_ribbon(aes(x = x, ymin = lower.ts - fit, ymax = upper.ts - fit), alpha = .2) +
    theme_bw() + 
    xlab("") + ylab("") +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"))
}

# Control Q-Q plot (i.e. standard design without bands)
ctrl_lineup <- function(dframe) {
  require(ggplot2)
  ggplot(aes(x = x, y = y), data = dframe) + 
    geom_smooth(aes(x, fit), colour="grey50", se=FALSE, method="loess") +
    geom_point() + 
    theme_bw() + 
    facet_wrap(~.sample, ncol=5) +
    labs(x = "", y = "") +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank(), plot.margin=grid::unit(c(0, 0, 0, 0), "cm"))
  
}



########################
# Simulation functions #
########################

qq_plot_info <- function(x, conf = .95, line = FALSE) {
  dh_bands <- sim_env(x, conf = conf, line = line)
  colnames(dh_bands)[4:5] <- paste(colnames(dh_bands)[4:5], "dh", sep = ".")
  ts_bands <- QQ.cb(x, alpha = 1 - conf, plot = FALSE)
  colnames(ts_bands)[3:4] <- paste(colnames(ts_bands)[3:4], "ts", sep = ".")
  RES <- dplyr::inner_join(x = dh_bands, y = ts_bands, by = c("x", "y"))
  return(RES)
}

sim_lineup <- function(n, nplots, mean = 0, sd = 1, conf = 0.95) {
  sims <- replicate(nplots, rnorm(n, mean = mean, sd = sd))
  sims <- as.data.frame(sims)
  
  qq_info <- lapply(sims, qq_plot_info, conf = conf, line = FALSE)
  names(qq_info) <- 1:length(qq_info)
  
  RES <- plyr::ldply(qq_info, function(df) df, .id = ".sample")
}

#######################
# Example usage       #
#######################

# NOTE: The above functions only work for standardized data, 
#       that is data with mean = 0 and sd = 1. So if your data
#       are in the vector x, then run scale(x) before using the
#       simulation function. Some day I will write it to be more
#       general... but probably not this summer.
# Simulate 20 q-q plots with 20 points from N(0, 1)
sim_data <- sim_lineup(n = 20, nplots = 20)
ctrl_lineup(sim_data)
std_lineup(sim_data)
std_ts_lineup(sim_data)
rot_lineup(sim_data)
rot_ts_lineup(sim_data)
rot2_lineup(sim_data)
rot2_ts_lineup(sim_data)
