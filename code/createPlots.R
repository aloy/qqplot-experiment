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

files <- dir("lineups/data")

for (fname in files) {
  dframe <- read.csv(sprintf("lineups/data/%s", fname))
  
  for (dsnfun in c("std_lineup", "rot_lineup", "rot2_lineup", "ctrl_lineup")) {
    fit2 <- dframe
    b <- (HLMdiag:::qqlineInfo(dframe$x[dframe$.n==20]))[2] # need motivation for that
    fit2$x[fit2$.n==20] <- dframe$x[dframe$.n==20]/b # change variance to 1, compare then
    idx <- grep("naive1.qq", names(fit2))
    fit2 <- ddply(fit2[,-idx], .(.n), transform, naive1.qq=qqnorm(x, plot.it=FALSE))
    idx <- grep("naive1.env", names(fit2))
    fit2 <- ddply(fit2[,-idx], .(.n), transform, naive1.env=sim_env(x))
    
    res <- parse_filename(gsub(".csv", "", fname))
    # outer panel  
    fit2$.sample <- fit2$.sample_outer
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("lineups/pdfs/%s",gsub(".csv","-outer.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("lineups/images%s", tmpfile), 
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
    file="lineups/picture-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("lineups/picture-details.csv"), append=TRUE)
    
    # inner panel  
    fit2$.sample <- fit2$.sample_inner
    eval(as.symbol(dsnfun))(fit2)    
    ggsave(file=sprintf("lineups/pdfs/%s",gsub(".csv","-inner.pdf", fname)))
    tmpfile <- sprintf("%s.svg",tempfile(tmpdir=""))
    lineup(dsnfun, fit2, filename= sprintf("lineups/images%s", tmpfile), 
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
    file="lineups/picture-details.csv", row.names=FALSE, sep=",",
    col.names=!file.exists("lineups/picture-details.csv"), append=TRUE)
  }


  
  
}
