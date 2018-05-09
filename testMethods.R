testMethods <- function(object, data, dataList, formula = "Y~1", print=FALSE, family=Gaussian()){
  
  cat(object$formulaFDboost, "\n")
  
  p <- length(object$baselearner)
  
  ### get coefficients
  cat("get coefficients", "\n")
  hatcoef <- try(coef(object))
  
  ### try the plot-functions
  cat("plot functions", "\n")
  par(mfrow=c(1,2))
  try(plotPredicted(object, ask = FALSE, lwdPred=2, subset=1:5))
  try(plotResiduals(object, ask = FALSE, subset=1:5))
  
  cat("plot coefficients", "\n")
  par(mfrow=c(3,p))
  try(plot(object, ask=FALSE, pers=TRUE))
  try(plot(object, ask=FALSE, pers=FALSE))
  try(plot(object, ask=FALSE, raw=TRUE))
  par(ask=FALSE)  
  
  cat("predict", "\n")
  #cat("predict effectwise toFDboost=TRUE", "\n")
  test0 <- predict(object)
  test1 <- predict(object, newdata=data)
  if(print){
    print(str(test0))
    print(str(test1))
  }
  for(i in 0:p){
    test0 <- try(predict(object, which=i))
    test1 <- try(predict(object, which=i, newdata=data))
    if(print){
      print(str(test0))
      print(str(test1))
    }
    if(any(test0 - test1 > 10^-6)){
      cat("ERROR: prediction is not equal for effet", i, "\n")
    }
  }
  test0 <- predict(object, which=1:p)
  test1 <- predict(object, newdata=data, which=1:p)
  if(print){
    print(str(test0))
    print(str(test1))
  }
  
  #cat("predict effectwise, toFDboost=FALSE", "\n")
  test0 <- predict(object, toFDboost=FALSE)
  test1 <- predict(object, newdata=data, toFDboost=FALSE)
  if(print){
    print(str(test0))
    print(str(test1))
  }
  for(i in 0:p){
    test0 <- try(predict(object, which=i, toFDboost=FALSE))
    test1 <- try(predict(object, which=i, newdata=data, toFDboost=FALSE))
    if(print){
      print(str(test0)) 
      print(str(test1)) 
    }
    if(any(test0 - test1 > 10^-6)){
      cat("ERROR: prediction is not equal for effet", i, "\n")
    }
  }
  test0 <- predict(object, which=1:p, toFDboost=FALSE)
  test1 <- predict(object, newdata=data, which=1:p, toFDboost=FALSE)
  if(print){
    print(str(test0))
    print(str(test1))
  }
  if(any(test0 - test1 > 10^-6)){
    cat("ERROR: prediction is not equal for effet 1:p", "\n")
  }
  
  
  for(i in 0:p){
    test0 <- try(predict(object, which=i, toFDboost=TRUE))
    test1 <- try(predict(object, which=i, newdata=data, toFDboost=TRUE))
    if(print){
      print(str(test0)) 
      print(str(test1)) 
    }
    if( any(test0 - test1 > 10^-6) ){
      cat("ERROR: prediction is not equal for effet", i, "\n")
    }
  }
  
  ### try the update-function
  cat("update function", "\n")
  #browser()
  newdat <- object$data
  resp <- object$response
  if(!is.factor(resp)) resp[1] <- 10
  nameResp <- as.character(as.formula(object$formulaFDboost)[[2]])
  newdat$resp <- if(!is.null(object$ydim)) I(matrix(resp, ncol=object$ydim[2])) else resp
  names(newdat)[names(newdat)=="resp"] <- nameResp
  newdat$Ynew <- newdat$Y
  # !!!!!!!!!!!!!! does not work in general when X1 / X1h are named differently !!!!!!!!!!!!!!
  # if(grepl("bhistx",object$formulaFDboost)){ 
  #   
  #   newdat$X1 <- I(attr(newdat$X1h, "x"))
  #   newdat$sTemp <- attr(newdat$X1h, "argvals")
  #   names(newdat)[names(newdat)=="sTemp"] <- attr(newdat$X1h, "argvalsLab")
  #   
  # }
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  newdat$timeHere <- object$yind 
  names(newdat)[names(newdat)=="timeHere"] <- attr(object$yind,"nameyind")
  newWei <- object$`(weights)`
  newWei[1:3] <- 0
  frmla <- paste0("Y ~ ", gsub("1 + ", "", as.character(as.formula(object$formulaFDboost))[[3]]))
  objectOrg <- object
  tryCatch(invisible(update(objectOrg, data=newdat)), 
           error=function(e)print(paste0("test data: ",e)),
           finally=print("test data END"))
  tryCatch(invisible(update(objectOrg, formula=update.formula(as.formula(objectOrg$formulaFDboost),"Ynew ~ ."),
                                    data=newdat)), 
           error=function(e)print(paste0("test formula+data: ",e)),
           finally=print("test formula+data END"))
  tryCatch(invisible(update(objectOrg, formula=update.formula(as.formula(objectOrg$formulaFDboost), 
                                                                   frmla)
  )), 
  error=function(e)print(paste0("test formula: ",e)),
  finally=print("test formula END"))
  tryCatch(invisible(update(objectOrg, weights=newWei)), 
           error=function(e)print(paste0("test weights: ",e)),
           finally=print("test weights END"))
  tryCatch(invisible(update(objectOrg, control=boost_control(mstop=3,nu=1))), 
           error=function(e)print(paste0("test control: ",e)),
           finally=print("test control END"))
  tryCatch(invisible(update(objectOrg, offset = "scalar")), error=function(e)print(paste0("test offset: ",e)),
           finally=print("test offset END"))
  
  
  ## try cvrisk / validateFDboost
  object <- object[40]
  
  cat("cvrsik()", "\n")
  cvm1 <- cvrisk(object, grid = 1:40, 
                          folds = cvLong(id = object$id, weights = model.weights(object), B=3))
  try(plot(cvm1))
  
  cat("validateFDboost()", "\n")
  val1 <- try(validateFDboost(object, grid = 1:40, 
                          folds = cv(rep(1, length(unique(object$id))), type = "bootstrap", B=3)))
  
  if(any(class(val1) != "try-error")){
    cat("plot validateFDboost", "\n")
    try(plot(val1))
    try(plotPredCoef(val1, terms=FALSE, which=length(object$basemodel), ask = FALSE))
    try(plotPredCoef(val1, terms=FALSE, pers=FALSE, which=length(object$basemodel), ask = FALSE))
  }
  
  # if( ! "FDboostScalar" %in% class(object) ){ ## applyFolds() suited to functional response 
    cat("applyFolds()", "\n")
    appl1 <- try(applyFolds(object, grid = 1:40, 
                            folds = cv(rep(1, length(unique(object$id))), type = "bootstrap", B=3)))
    if(any(class(appl1) != "try-error")){ 
      try(plot(appl1, ask = FALSE))
    }
  # }

    cat("bootstrapCI()", "\n")
    bsci1 <- try(bootstrapCI(object, B_inner = 3, B_outer = 5))
    if(any(class(bsci1) != "try-error")){ 
      try(plot(bsci1, ask = FALSE))
      # old function:
      # par(mfrow=c(2,2))
      # try(plot_bootstrapCI(bsci1, commonRange = FALSE, ask = FALSE))
      # try(plot_bootstrapCI(bsci1, commonRange = FALSE, pers = FALSE, ask = FALSE))
    }

  #return(list(pred0=test0, pred1=test1))
  
}