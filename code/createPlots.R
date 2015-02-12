#-------------------------------------------------------------------------------
# Script generating lineups of Q-Q plots: standard and detrended versions
#
# Heike Hofmann, Adam Loy
# February 2015
#
# References: 
#   Aldor-Noiman, S., Brown, L.D., Buja, A., Rolke, W., and Stine, R.A. (2013), 
#   “The Power to See: A New Graphical Test of Normality,” The American Statistician, 
#   67, 249–260.
#-------------------------------------------------------------------------------


setwd("~/Dropbox/lineups-nf/turk17/")
library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(mlmRev)   # for the exam data
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(gridSVG)

# required to TS test bands
library(mvtnorm)    # for the multivariate normal distribution
library(robustbase) # for robust estimates for the mean and sd


scriptURL = "http://www.hofroe.net/examples/lineup/qqaction.js"

make_interactive <- function(filename, script, toggle="toggle", background="#e5e5e5") {
# use #ffffff for white background
  require(gridSVG)
  grobs <- grid.ls()
  
  idx <- grep("panel-", grobs$name)
  for (i in idx) { 
    grid.garnish(grobs$name[i],
                 onmouseover=paste("frame('",grobs$name[i+2], ".1')", sep=""),
                 onmouseout=paste("deframe('",grobs$name[i+2], ".1')", sep=""), 
                 onmousedown=sprintf("%shigh(evt, '%s.1', '%s')", toggle, grobs$name[i+2], background)
                 )
  }

  # use script on server to get locally executable javascript code
  # or use inline option
  grid.script(filename=script)
  grid.export(filename, uniqueNames=FALSE, exportJS="inline", exportCoords="inline", exportMappings="inline")
}

std_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  print(ggplot(aes(x=naive1.qq.x, y=naive1.qq.y), data = dframe) + 
          #  geom_smooth(aes(naive1.qq.x, naive1.env.fit.value), colour="grey50", se=FALSE, method="loess") +
          geom_abline(colour="grey50") +
          facet_wrap(~.sample, ncol=5) +
          geom_point() + 
          theme_bw() + xlab("") + ylab("") +
          geom_ribbon(aes(x = naive1.env.fit.value, ymin = naive1.env.lower, ymax = naive1.env.upper),alpha = .2)+
          theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                axis.ticks = element_blank(), plot.margin=unit(c(0,0,0,0), "cm"))
  )  
}

std_ts_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  print(ggplot(aes(x=ts.qq.x, y=ts.qq.y), data = dframe) + 
          #  geom_smooth(aes(naive1.qq.x, naive1.env.fit.value), colour="grey50", se=FALSE, method="loess") +
          geom_abline(colour="grey50") +
          facet_wrap(~.sample, ncol=5) +
          geom_point() + 
          theme_bw() + xlab("") + ylab("") +
          geom_ribbon(aes(x = ts.qq.x, ymin = ts.qq.lower, ymax = ts.qq.upper),alpha = .2)+
          theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                axis.ticks = element_blank(), plot.margin=unit(c(0,0,0,0), "cm"))
  )  
}

rot2_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  print(  ggplot(aes(x=naive1.qq.x, y=naive1.qq.y-naive1.env.fit.value), data = dframe) + 
            facet_wrap(~.sample, ncol=5) + 
            geom_hline(yintercept=0, colour="grey30")+
            geom_point() + 
            geom_ribbon(aes(x = naive1.qq.x, 
                            ymin = naive1.env.lower-naive1.env.fit.value, ymax = naive1.env.upper-naive1.env.fit.value),alpha = .2)+
            theme_bw() + xlab("") + ylab("") +
            ylim(range(dframe$naive1.qq.y)) +
            theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                  axis.ticks = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")))
}

rot2_ts_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  print(  ggplot(aes(x=ts.qq.x, y=ts.qq.y-ts.qq.x), data = dframe) + 
            facet_wrap(~.sample, ncol=5) + 
            geom_hline(yintercept=0, colour="grey30")+
            geom_point() + 
            geom_ribbon(aes(x = ts.qq.x, 
                            ymin = ts.qq.lower-ts.qq.x, ymax = ts.qq.upper-ts.qq.x),alpha = .2)+
            theme_bw() + xlab("") + ylab("") +
            ylim(range(dframe$ts.qq.y)) +
            theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                  axis.ticks = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")))
}

rot_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  print(  ggplot(aes(x=naive1.qq.x, y=naive1.qq.y-naive1.env.fit.value), data = dframe) + 
            facet_wrap(~.sample, ncol=5) + 
            geom_hline(yintercept=0, colour="grey30")+
            geom_point() + 
            geom_ribbon(aes(x = naive1.qq.x, 
                            ymin = naive1.env.lower-naive1.env.fit.value, ymax = naive1.env.upper-naive1.env.fit.value),alpha = .2)+
            theme_bw() + xlab("") + ylab("") +
#            xlab("Normal Quantiles") + 
#            ylab("Sample Quantiles") +
            theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                  axis.ticks = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")))
}

rot_ts_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  print(  ggplot(aes(x=ts.qq.x, y=ts.qq.y-ts.qq.x), data = dframe) + 
            facet_wrap(~.sample, ncol=5) + 
            geom_hline(yintercept=0, colour="grey30")+
            geom_point() + 
            geom_ribbon(aes(x = ts.qq.x, 
                            ymin = ts.qq.lower-ts.qq.x, ymax = ts.qq.upper-ts.qq.x),alpha = .2)+
            theme_bw() + xlab("") + ylab("") +
            theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                  axis.ticks = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")))
}

ctrl_lineup <- function(dframe) {
  require(ggplot2)
  require(grid)
  print(ggplot(aes(x=naive1.qq.x, y=naive1.qq.y), data = dframe) + 
          geom_smooth(aes(naive1.qq.x, naive1.env.fit.value), colour="grey50", se=FALSE, method="loess") +
          geom_point() + 
          facet_wrap(~.sample, ncol=5) +
          theme_bw() + xlab("") + ylab("") +
          theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                axis.ticks = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")) 
  )
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
  
  q.prob<-qnorm((1:n)/(n+1))
  z.sample<-(x.sample-mu)/sigma
  if(plot==TRUE){
    plot(q.prob,norm.upper,type="l",col="red",ylim=c(-3,3),xlim=c(-3,3),ylab="Sample Quantile",xlab="Sample Quantile")
    lines(q.prob,norm.lower,col="red")
    points(q.prob,sort(z.sample),pch=19,cex=0.6)
  }
  
  return(data.frame(x=q.prob, y=sort(z.sample), lower=norm.lower, upper=norm.upper))
}


sim_env <- function(x, conf = .95, line=FALSE){
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
  return(data.frame(lower, upper, fit.value, idx=1:length(x)))
  
}

parse_filename <- function(fname) {
	pars <- strsplit(fname, "-")[[1]]
	names(pars) <- c("data", "null", "rep",  "i", "n", "df", "innerPos", "outerPos")
	pars
}


lineup <- function(method, dframe, filename, script, toggle="toggle") {
  z = as.list(match.call()[-1])
	
  eval(call(eval(z$method), dframe))
  require(gridSVG)
  grobs <- grid.ls()
  
  idx <- grep("panel-", grobs$name)
  for (i in idx) { 
    grid.garnish(grobs$name[i],
                 onmouseover=paste("frame('",grobs$name[i+2], ".1')", sep=""),
                 onmouseout=paste("deframe('",grobs$name[i+2], ".1')", sep=""), 
                 onmousedown=paste(sprintf("%shigh(evt, '", toggle),grobs$name[i+2], ".1')", sep=""))
  }

  # use script on server to get locally executable javascript code
  # or use inline option
  grid.script(filename=script)
  grid.export(filename, uniqueNames=FALSE, exportJS="inline", exportCoords="inline", exportMappings="inline")
}
##########################
# for N(0,s)
files <- dir("lineups-N(0,s)/data")

for (fname in files) {
  dframe <- read.csv(sprintf("lineups-N(0,s)/data/%s", fname))
  
  fit2 <- dframe
  b <- (HLMdiag:::qqlineInfo(dframe$x[dframe$.n==20]))[2] # slope of the real data determines slope of all
  fit2$x[fit2$.n==20] <- dframe$x[dframe$.n==20]/b # change variance to 1, compare then
  idx <- grep("naive1.qq", names(fit2))
  fit2 <- ddply(fit2[,-idx], .(.n), transform, 
                naive1.qq=qqnorm(x, plot.it=FALSE))
  idx <- grep("naive1.env", names(fit2))
  fit2 <- ddply(fit2[,-idx], .(.n), transform, 
                naive1.env=sim_env(x),
                ts.qq=QQ.cb(x, plot=FALSE))
  res <- parse_filename(gsub(".csv", "", fname))

  for (dsnfun in c( "std_ts_lineup", "rot_ts_lineup", "rot2_lineup", "rot2_ts_lineup")) {
    # outer panel  
    fit2$.sample <- fit2$.sample_outer
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("lineups-N(0,s)/pdfs/%s",gsub(".csv","-outer.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("lineups-N(0,s)/images%s", tmpfile), 
           toggle="toggle", script=scriptURL)
    test_param <- sprintf("turk17-%s-Outer-Multiple-rep-%s-null-%s", dsnfun, res["rep"], res["null"])
    param_value <- sprintf("%s-%s", res["n"], res["df"])

    write.table(data.frame(
      sample_size=res["n"], 
      test_param=test_param,
      param_value=param_value,
      p_value=NA,
      obs_plot_location=res["outerPos"], 
      pic_name=tmpfile,
      experiment="turk17",
      difficulty=fname,
      data_name=fname
    ), 
    file="lineups-N(0,s)/picture-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("lineups-N(0,s)/picture-details.csv"), append=TRUE)
    
    # inner panel  
    fit2$.sample <- fit2$.sample_inner
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("lineups-N(0,s)/pdfs/%s",gsub(".csv","-inner.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("lineups-N(0,s)/images%s", tmpfile), 
           toggle="toggle", script=scriptURL)
    test_param <- sprintf("turk17-%s-Inner-Multiple-rep-%s-null-%s", dsnfun, res["rep"], res["null"])
    param_value <- sprintf("%s-%s", res["n"], res["df"])
    
    write.table(data.frame(
      sample_size=res["n"], 
      test_param=test_param,
      param_value=param_value,
      p_value=NA,
      obs_plot_location=res["innerPos"], 
      pic_name=tmpfile,
      experiment="turk17",
      difficulty=fname,
      data_name=fname
    ), 
    file="lineups-N(0,s)/picture-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("lineups-N(0,s)/picture-details.csv"), append=TRUE)
  }


  
  
}
#####################################
files <- dir("lineups-N(0,1)/data")

for (fname in files) {
  dframe <- read.csv(sprintf("lineups-N(0,1)/data/%s", fname))
  
  fit2 <- dframe
#  b <- (HLMdiag:::qqlineInfo(dframe$x[dframe$.n==20]))[2] # slope of the real data determines slope of all
#  fit2$x[fit2$.n==20] <- dframe$x[dframe$.n==20]/b # change variance to 1, compare then
  idx <- grep("naive1.qq", names(fit2))
  fit2 <- ddply(fit2[,-idx], .(.n), transform, 
                naive1.qq=qqnorm(x, plot.it=FALSE))
  idx <- grep("naive1.env", names(fit2))
  fit2 <- ddply(fit2[,-idx], .(.n), transform, 
                naive1.env=sim_env(x),
                ts.qq=QQ.cb(x, plot=FALSE))
  res <- parse_filename(gsub(".csv", "", fname))
  
for (dsnfun in c( "std_ts_lineup", "rot_ts_lineup", "rot2_lineup", "rot2_ts_lineup")) {
  # outer panel  
    fit2$.sample <- fit2$.sample_outer
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("lineups-N(0,1)/pdfs/%s",gsub(".csv","-outer.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("lineups-N(0,1)/images%s", tmpfile), 
           toggle="toggle", script=scriptURL)
    test_param <- sprintf("turk17-%s-Outer-Multiple-rep-%s-null-%s", dsnfun, res["rep"], res["null"])
    param_value <- sprintf("%s-%s", res["n"], res["df"])
    
    write.table(data.frame(
      sample_size=res["n"], 
      test_param=test_param,
      param_value=param_value,
      p_value=NA,
      obs_plot_location=res["outerPos"], 
      pic_name=tmpfile,
      experiment="turk17",
      difficulty=fname,
      data_name=fname
    ), 
    file="lineups-N(0,1)/picture-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("lineups-N(0,1)/picture-details.csv"), append=TRUE)
    
    # inner panel  
    fit2$.sample <- fit2$.sample_inner
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("lineups-N(0,1)/pdfs/%s",gsub(".csv","-inner.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("lineups-N(0,1)/images%s", tmpfile), 
           toggle="toggle", script=scriptURL)
    test_param <- sprintf("turk17-%s-Inner-Multiple-rep-%s-null-%s", dsnfun, res["rep"], res["null"])
    param_value <- sprintf("%s-%s", res["n"], res["df"])
    
    write.table(data.frame(
      sample_size=res["n"], 
      test_param=test_param,
      param_value=param_value,
      p_value=NA,
      obs_plot_location=res["innerPos"], 
      pic_name=tmpfile,
      experiment="turk17",
      difficulty=fname,
      data_name=fname
    ), 
    file="lineups-N(0,1)/picture-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("lineups-N(0,1)/picture-details.csv"), append=TRUE)
  }
  
  
  
  
}

##########################
# trials

files <- dir("trials/data")

for (fname in files) {
  dframe <- read.csv(sprintf("trials/data/%s", fname))
  
  fit2 <- dframe
  b <- (HLMdiag:::qqlineInfo(dframe$x[dframe$.n==20]))[2] # slope of the real data determines slope of all
  fit2$x[fit2$.n==20] <- dframe$x[dframe$.n==20]/b # change variance to 1, compare then
  idx <- grep("naive1.qq", names(fit2))
  fit2 <- ddply(fit2[,-idx], .(.n), transform, 
                naive1.qq=qqnorm(x, plot.it=FALSE))
  idx <- grep("naive1.env", names(fit2))
  fit2 <- ddply(fit2[,-idx], .(.n), transform, 
                naive1.env=sim_env(x),
                ts.qq=QQ.cb(x, plot=FALSE))
  res <- parse_filename(gsub(".csv", "", fname))
  
  #  for (dsnfun in c("std_lineup", "std_ts_lineup", "rot_lineup", "rot_ts_lineup", "rot2_lineup", "rot2_ts_lineup", "ctrl_lineup")) {
  for (dsnfun in c( "std_ts_lineup", "rot_ts_lineup", "rot2_lineup", "rot2_ts_lineup")) {
    # outer panel  
    fit2$.sample <- fit2$.sample_outer
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("trials/pdfs/%s",gsub(".csv","-outer.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("trials/images%s", tmpfile), 
           toggle="toggle", script=scriptURL)
    test_param <- sprintf("turk17-%s-Outer-Multiple-rep-%s-null-%s", dsnfun, res["rep"], res["null"])
    param_value <- sprintf("%s-%s", res["n"], res["df"])
    
    write.table(data.frame(
      sample_size=res["n"], 
      test_param=test_param,
      param_value=param_value,
      p_value=NA,
      obs_plot_location=res["outerPos"], 
      pic_name=tmpfile,
      experiment="turk17",
      difficulty=fname,
      data_name=fname
    ), 
    file="trials/picture-trials-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("trials/picture-trials-details.csv"), append=TRUE)
    
    # inner panel  
    fit2$.sample <- fit2$.sample_inner
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("trials/pdfs/%s",gsub(".csv","-inner.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("trials/images%s", tmpfile), 
           toggle="toggle", script=scriptURL)
    test_param <- sprintf("turk17-%s-Inner-Multiple-rep-%s-null-%s", dsnfun, res["rep"], res["null"])
    param_value <- sprintf("%s-%s", res["n"], res["df"])
    
    write.table(data.frame(
      sample_size=res["n"], 
      test_param=test_param,
      param_value=param_value,
      p_value=NA,
      obs_plot_location=res["innerPos"], 
      pic_name=tmpfile,
      experiment="turk17",
      difficulty=fname,
      data_name=fname
    ), 
    file="trials/picture-trials-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("trials/picture-trials-details.csv"), append=TRUE)
  }
  
  
  
  
}

