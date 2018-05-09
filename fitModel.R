fitModel <- function(formula, scenario, n = 100, nxgrid = 35, nygrid = 60, SNR = 10,
                     propmissing = 0, limits = NULL, family = Gaussian(), ...){

  dots <- list(...)
  
  #   data1 <- pffrSim(scenario = c( "int", "ff"), n = 100, nxgrid = 40, nygrid = 60, SNR = 10,
  #                    propmissing = 0, limits = NULL)
  
  data1 <- pffrSim(scenario = scenario, n = n, nxgrid = nxgrid, nygrid = nygrid, SNR = SNR,
                   propmissing = propmissing, limits = limits)
  
  if(family@name != "Squared Error (Regression)"){
    if(family@name %in% c("Negative Binomial Likelihood", "Adaboost Exponential Error")){
      data1$Y0 <- data1$Y
      data1$Y <- I(matrix("one", ncol = ncol(data1$Y), nrow = nrow(data1$Y) ))
      data1$Y[ data1$Y0 > mean(data1$Y0) ] <- "two"
    }
    
    if(family@name %in% c("Negative Gamma Likelihood")){
      data1$Y <- I(abs(data1$Y))
    }
    
    if(family@name %in% c("Poisson Likelihood")){
      #data1$Y <- I(round(abs(data1$Y)))
      data1$Y <- apply(exp(data1$Y/2), 2, rpois, n = nrow(data1$Y) )
    }
  }
  
  if(propmissing==0){
    dataList <- as.list(data1)
    dataList$tvals <- attr(data1, "yindex")
    
    X1h <- with(dataList, hmatrix(time= rep(tvals, each=n), id=rep(1:n, nygrid), 
                                  x=X1, argvals=attr(data1, "xindex"), 
                                  timeLab="tvals", idLab="wideIndex", xLab="myX", argvalsLab="svals"))
    dataList$X1h <- I(X1h)           
    
  }else{
    dataList <- as.list(data1$data)
    dataList$Y <- data1$ydata$.value[order(data1$ydata$.index)]
    dataList$idy <- data1$ydata$.obs[order(data1$ydata$.index)]
    dataList$tvals <- data1$ydata$.index[order(data1$ydata$.index)]
    attr(dataList, "truth") <- attr(data1, "truth")
    attr(dataList, "xindex") <- attr(data1, "xindex")
    attr(dataList, "yindex") <- attr(data1, "yindex")
    
    X1h <- with(dataList, hmatrix(time= tvals, id=idy, 
                                  x=X1, argvals=attr(data1, "xindex"), 
                                  timeLab="tvals", idLab="wideIndex", xLab="myX", argvalsLab="svals"))
    dataList$X1h <- I(X1h)
    
  } 
  
  dataList$svals <- attr(data1, "xindex")
  
  ## center the data when necessary
  dataList$X1 <- scale(dataList$X1, scale=FALSE)
  dataList$xlin <- scale(dataList$xlin, scale=FALSE)
  
  dataList$zlong <- factor(gl(n=2, k=n/2, length=n*nygrid*(1-propmissing)), levels=1:3)  ## add a factor variable
  dataList$z <- factor(gl(n=2, k=n/2, length=n), levels=1:3)
  dataList$zlong2 <- gl(n=2, k=1, length=n*nygrid*(1-propmissing))  ## add a factor variable
  dataList$z2 <- gl(n=2, k=1, length=n)
  ## formula=Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5) %X2% bolsc(zlong)
  ## formula=Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5) %X2% bolsc(zlong) %X2% bolsc(zlong2)
  
  ## formula=Y ~ 1 + bhistx(x = X1h, df = 5)
  ## formula=Y ~ 1 + bhistx(x = X1h, df = 5) %X% bolsc(zlong)
  
  if(propmissing==0){ 
    
    if(nygrid > 1){ ## functional response 
      
      ## fit the model in wide format
      mod <- FDboost(formula, timeformula=~bbs(tvals, knots=10, df=3), data=dataList, family=family, ...)  ##, ...)
      ## mod <- FDboost(formula, timeformula=~bbs(tvals, knots=10), data=dataList) 
      #if(any(grepl("bhist", formula))){
      #  ## compare hist effect to truth
      #  funplot(mod$yind, attr(dataList, "truth")$etaTerms$X1[1:5,], lwd=2)
      #  funplot(mod$yind, predict(mod, which=2)[1:5,], add=TRUE) 
      #} 
      ## pred <- predict(mod)
      ## pred <- predict(mod, newdata=dataList[c("X1","svals","tvals", "zlong")], which=2)
      ## pred <- predict(mod, newdata=list(X1h=dataList[["X1h"]], zlong=dataList[["zlong"]], tvals=attr(data1, "yindex")), which=2)
      ##pred <- predict(mod, newdata=list(X1h=dataList[["X1h"]], zlong=dataList[["zlong"]], tvals=attr(data1, "yindex")), which=2)
      
      ## pred <- predict(mod, newdata=dataList[c("X1","svals","tvals", "zlong", "zlong2")], which=2)
      ##pred <- predict(mod, newdata=list(X1=I(dataList$X1[1:20,1:20]), 
      ##                                  tvals=dataList$tvals[1:20], svals=dataList$svals[1:20]), which=2)
      ##pred <- predict(mod, newdata=data.frame(X1=I(dataList$X1[1:20,1:20]), 
      ##                                  tvals=dataList$tvals[1:20], svals=dataList$svals[1:20], hello=1:20), which=2)
      
      ## preda <- predict(mod, which=2)
      ## cf <- coef(mod, which=2)
      ## plot(mod, which=2)
      
    }else{ ## scalar response as nygrid == 1
      
      mod <- FDboost(formula, timeformula=NULL, data=dataList, family=family, ...)
      
    }

  }else{   
    ## fit the model in long format
    mod <- FDboost(formula, 
                   timeformula=~bbs(tvals, knots=10, df=3), id=~idy, data=dataList, family=family, ...)
    ## pred <- predict(mod, newdata=dataList[c("X1","svals","tvals", "zlong", "zlong2", "idy")], which=2)
    ## pred <- predict(mod, newdata=list(X1h=dataList[["X1h"]], zlong=dataList$zlong, tvals=dataList$tvals), which=2)
    
    if(FALSE){
      myinS <- "linear"
      mod <- FDboost(Y ~ 1 + bsignal(x = X1, s = svals, df = 5, inS = myinS), 
                     timeformula=~bbs(tvals, knots=10, df=3), id=~idy, data=dataList)
    }
    
    #if(any(grepl("bhist", formula))){
    #  ## compare hist effect to truth
    #  funplot(attr(dataList, "yindex"), attr(dataList, "truth")$etaTerms$X1[1:5,], lwd=2)
    #  funplot(mod$yind[mod$id %in% 1:5], predict(mod, which=2)[mod$id %in% 1:5], 
    #          mod$id[mod$id %in% 1:5], add=TRUE)
    #} 
  }
  
  return(list(mod=mod, dat=dataList))
  
}