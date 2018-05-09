
# library(FDboost)
library(splines)
## library(refund)

# setwd("Y:/transfer/SarahDavid/test_FDboost")


#devtools::install_github("refunders/refund", ref="devel")
#library(refundDevel)
library(refund)
### source("../functions/update.FDboost.R")
source("fitModel.R")
source("testMethods.R")

# ## generate the data 
# data2 <- pffrSim(scenario = "all", n = 100, nxgrid = 40, nygrid = 60, SNR = 10,
#         propmissing = 0, limits = NULL)
# 
# ## center the data where necessary
# data2$X1 <- scale(data2$X1, scale=FALSE)
# data2$xlin <- data2$xlin - mean(data2$xlin)
# 
# dataList <- as.list(data2)
# dataList$tvals <- attr(data2, "yindex")
# dataList$svals <- attr(data2, "xindex")
# 
# ## fit the model
# m2 <- FDboost(Y ~ bsignal(x=X1, s=svals, df=5) + #linear function-on-function
#              bolsc(xlin, df=5)  +  #varying coefficient term
#              c(bbs(xte1, xte2, df=5)) + #bivariate smooth term in xte1 & xte2, const. over Y-index
#              bbsc(xsmoo, df=5, knots=10) + #smooth effect of xsmoo varying over Y-index
#              c(bols(xconst, df=5)), # linear effect of xconst constant over Y-index
#            timeformula=~bbs(tvals, knots=10),
#            data=dataList)
# 
# extract(m2, "df")


########################################## 

### simulate the effects one by one, always include an intercept

if(FALSE){
  formula=Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5)
  ## formula= Y ~ 1 + bsignal(x=X1, s=svals, df=5)
  scenario=c("int", "ff") 
  propmissing=0
  limits = function(s,t){ s<=t }
  n = 100; nxgrid = 35; nygrid = 60; SNR = 10
}

if(FALSE){
  formula=Y ~ 1 + bsignal(x=X1, s=svals, df=5)
  scenario=c("int", "ff") 
  propmissing=0
  limits = NULL
  n = 100; nxgrid = 35; nygrid = 60; SNR = 10
}

if(FALSE){
  
  s <- seq(0, 1, length = 20)
  t <- seq(0, 1, length = 20)
  
  test1 <- function(s, t) {
    s * cos(pi * abs(s - t)) - 0.19
  }
  
  beta1.st <- outer(s, t, test1)
  #   if (!is.null(limits)) {
  #     range <- outer(s, t, limits)
  #     beta1.st <- beta1.st * range
  #   }
  
  par(mfrow=c(1,2))
  persp(beta1.st, theta=30, phi=30, ticktype="detailed")
  plot(mod, which=2, pers=T)
  
}

options(warn = 1)

############# pure intercept model 
m0 <- fitModel(formula = Y ~ 1, scenario=c("int"))
testMethods(object=m0$mod, data=m0$dat, dataList = m0$dat, formula=Y ~ 1)

m0i <- fitModel(formula=Y ~ 1, scenario=c("int"), propmissing=0.3)
testMethods(object=m0i$mod, data=m0i$data, dataList=m0i$dat, formula=Y ~ 1)


############# linear function-on-function,  \int x(s)beta(s,t)ds
m1 <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"))
testMethods(object=m1$mod, data=m1$dat, dataList = m1$dat, formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3))
#pred0 <- predict(m1$mod)
#pred <- predict(m1$mod, newdata=m1$dat)
#preddat <- predict(m1$mod, newdata=data.frame(X1=m1$dat$X1[1:35,], svals=m1$dat$svals, tvals=m1$dat$tvals[1:35]) )
## all(preddat==pred[1:35, 1:35])

m1i <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), propmissing=0.3)
## gives error in prediction, as data for model fit cannot be directly used for prediciton,
## as id is used in the model fit but NOT in prediction
# testMethods(m1i$mod, data=m1i$dat)
## constuct the according data for the prediction 
newd <- list(X1=m1i$dat$X1[m1i$dat$idy,], svals=m1i$dat$svals, tvals=m1i$dat$tvals)
testMethods(object=m1i$mod, data=newd, dataList=m1i$dat, formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3))
pred0 <- predict(m1i$mod)
pred <- predict(m1i$mod, newdata=newd)
pred <- predict(m1i$mod, newdata=newd, which=2)
# pred <- predict(m1i$mod, newdata=m1i$dat, which=2)
preddat <- predict(m1i$mod, newdata=data.frame(X1=m1i$dat$X1[1:35,], svals=m1i$dat$svals, tvals=m1i$dat$tvals[1:35]) )
## cf <- coef(m1i$mod)
## plot(m1i$mod, pers=TRUE)

### test with further arguments
KNOTS <- 5
mydf <- 3
m1ia <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=mydf, knots=KNOTS), scenario=c( "int", "ff"), 
                propmissing=0.3, numInt="Riemann")
extract(m1ia$mod, "df")
dim(extract(m1ia$mod, "design")[[2]])
dim(extract(m1i$mod, "design")[[1]])

newda <- list(X1=m1ia$dat$X1[m1ia$dat$idy,], svals=m1ia$dat$svals, tvals=m1ia$dat$tvals)
testMethods(object=m1ia$mod, data=newda, dataList=m1ia$dat, formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3, knots=KNOTS))


############# historical linear function-on-function,  \int_0^t x(s)beta(s,t)ds
m2 <- fitModel(formula = Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5), 
               scenario=c("int", "ff"), limits = function(s,t){ s<=t }, nxgrid=39)
testMethods(object=m2$mod, data=m2$dat, dataList = m2$dat, formula=Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5))

# ## compare hist effect to truth
# funplot(m2$mod$yind, attr(m2$dat, "truth")$etaTerms$X1[1:5,])
# funplot(m2$mod$yind, predict(m2$mod, which=2)[1:5,], add=TRUE, lwd=2)

### test with further arguments, in particual with limits-function
KNOTS = 25
m2a <- fitModel(formula= Y ~ 1 + bhist(x=X1, s=svals, time=tvals, knots = KNOTS, df=5, limits = function(xx, yy) xx <= yy ), 
                scenario=c( "int", "ff"), numInt="Riemann")
testMethods(object=m2a$mod, data=m2a$dat, dataList=m2a$dat, 
            formula= Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5, limits = function(xx, yy) xx <= yy ))
###


m2i <- fitModel(formula= Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5), 
                scenario=c( "int", "ff"), limits = function(s,t){ s<=t }, propmissing=0.3, nxgrid=39)
testMethods(object=m2i$mod, data=m2i$dat, dataList = m2i$dat, formula=Y ~ 1 + bhist(x=X1, s=svals, time=tvals, df=5))
## cf <- coef(m2i$mod)
newd <- list(X1=m2i$dat$X1[m2i$dat$idy,], svals=m2i$dat$svals, tvals=m2i$dat$tvals, idy=1:length(m2i$dat$tvals))
pred <- predict(m2i$mod, newdata=newd, which=2)
## predh <- predict(m2i$mod, newdata=m2i$dat, which=2)
## preddat <- predict(m2i$mod, newdata=data.frame(X1=m2i$dat$X1[1:39,], svals=m2i$dat$svals, tvals=m2i$dat$tvals[1:39]) )
## all(preddat==pred[1:39, 1:39])

# ## hist effect
# funplot(m2$mod$yind, attr(m2i$dat, "truth")$etaTerms$X1[1:5,], lwd=2)
# funplot(m2i$mod$yind[m2i$mod$id %in% 1:5], predict(m2i$mod, which=2)[m2i$mod$id %in% 1:5], 
#         m2i$mod$id[m2i$mod$id %in% 1:5], add=TRUE)



############# historical linear function-on-function,  \int_0^t x(s)beta(s,t)ds
m3 <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=5), 
               scenario=c("int", "ff"), limits = function(s,t){ s<=t }, nxgrid=39)
testMethods(object=m3$mod, data=m3$dat, dataList = m3$dat, formula=Y ~ 1 + bhistx(x=X1h, df=5))

test <- coef(m3$mod, which=2)

## further arguments
m3a <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=5, limits = function(xx,yy){ xx <= yy }), 
               scenario=c("int", "ff"), limits = function(s,t){ s <= t }, nxgrid=39)
testMethods(object=m3a$mod, data=m3a$dat, dataList = m3a$dat, formula=Y ~ 1 + bhistx(x=X1h, df=5))

# ## compare hist effect to truth
# funplot(m3$mod$yind, attr(m3$dat, "truth")$etaTerms$X1[1:5,])
# funplot(m3$mod$yind, predict(m3$mod, which=2)[1:5,], add=TRUE, lwd=2)

m3i <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=5), propmissing=0.3, 
               scenario=c("int", "ff"), limits = function(s,t){ s<=t }, nxgrid=39)
testMethods(object=m3i$mod, data=m3i$dat, dataList = m3i$dat, formula=Y ~ 1 + bhistx(x=X1h, df=5))
plot(m3i$mod, ask=FALSE)
test <- coef(m3i$mod, which=2)


############# historical linear function-on-function,  \int_0^t x(s)beta(s,t)ds
m4 <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=3) %X% bolsc(zlong, df=2), 
               scenario=c("int", "ff"), limits = function(s,t){ s<=t }, nxgrid=39)
test <- coef(m4$mod, which=2)
test$smterms[[1]]$main

testMethods(object=m4$mod, data=m4$dat, dataList = m4$dat, 
            formula=Y ~ 1 + bhistx(x=X1h, df=3) %X% bolsc(zlong, df=2))

test <- coef(m4$mod, which=2)

# ## compare hist effect to truth
# funplot(m3$mod$yind, attr(m3$dat, "truth")$etaTerms$X1[1:5,])
# funplot(m3$mod$yind, predict(m3$mod, which=2)[1:5,], add=TRUE, lwd=2)

m4i <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% bolsc(zlong, df=2), propmissing=0.3, 
                scenario=c("int", "ff"), limits = function(s,t){ s<=t }, nxgrid=39)
testMethods(object=m4i$mod, data=m4i$dat, dataList = m4i$dat, 
            formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% bolsc(zlong, df=2))
plot(m4i$mod)
test <- coef(m4i$mod, which=2)

m1$dat$index <- rep(1:nrow(m1$dat$X1), length(m1$dat$tvals))
myBlg <- with(m1$dat, bols(z, index = index) %Xc% 
                bols(z2, index = index))

m4ii <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% myBlg, 
                scenario=c("int", "ff"), nxgrid=39)
testMethods(object=m4ii$mod, data=m4ii$dat, dataList = m4ii$dat, 
            formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% myBlg)
plot(m4ii$mod)
test <- coef(m4ii$mod, which=2)

## use index in base-learner 
index <- m1$dat$index

m4iii <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% bolsc(z, index = index, df=2), 
                 scenario=c("int", "ff"), nxgrid=39)
testMethods(object=m4iii$mod, data=m4iii$dat, dataList = m4iii$dat, 
            formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% bolsc(z, index = index, df=2))
plot(m4iii$mod)
test <- coef(m4iii$mod, which=2)

m4iv <- fitModel(formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% bolsc(z, index = index, df=2), 
                 propmissing=0.00000001, 
                scenario=c("int", "ff"), limits = function(s,t){ s<=t }, nxgrid=39)
testMethods(object=m4iv$mod, data=m4iv$dat, dataList = m4iv$dat, 
            formula=Y ~ 1 + bhistx(x=X1h, df=5) %X% bolsc(z, index = index, df=2))
plot(m4iv$mod)
test <- coef(m4iv$mod, which=2)


############# model with smooth scalar term 
m5 <- fitModel(formula=Y ~ 1 + bbsc(xsmoo, df = 3), scenario=c( "int", "smoo"))
testMethods(object = m5$mod, data=m5$dat, dataList = m5$dat, formula=Y ~ 1 + bbsc(xsmoo, df = 3))

m5i <- fitModel(formula=Y ~ 1 + bbsc(xsmoo, df = 3), propmissing=0.3, scenario=c("int", "smoo"))
## constuct the according data for the prediction 
newd <- list(xsmoo=m5i$dat$xsmoo[m5i$dat$idy], tvals=m5i$dat$tvals)
testMethods(object=m5i$mod, data=newd, dataList = m5i$dat, 
            formula=Y ~ 1 + bbsc(xsmoo, df = 3))


## use c() notation
m5c <- fitModel(formula=Y ~ 1 + c(bbsc(xsmoo, df=2)), scenario=c( "int", "smoo"))
testMethods(m5c$mod, data=m5c$dat, dataList = m5c$dat, formula=Y ~ 1 + c(bbsc(xsmoo, df=2)))

## use c() notation with bolsc
m5c2 <- fitModel(formula=Y ~ 1 + c(bolsc(xsmoo, df=1)), scenario=c( "int", "smoo"))
testMethods(m5c2$mod, data=m5c2$dat, dataList = m5c2$dat, formula=Y ~ 1 + c(bolsc(xsmoo, df=1)))



############# model with concurrent term 
m6 <- fitModel(formula=Y ~ 1 + bconcurrent(x=X1, s=svals, time=tvals, df=3), 
               scenario=c("int", "ff"), limits = function(s,t){ s==t }, nxgrid=39, nygrid=39)
testMethods(object=m6$mod, data=m6$dat, dataList = m6$dat, formula=Y ~ 1 + bconcurrent(x=X1, s=svals, time=tvals, df=3))



############# models with scalar response 

## sof/ scalar response with functional effect 
ms1 <- fitModel(formula =  Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = Gaussian(), nygrid = 1)  ## scalar response as functional response with 1 observation 
testMethods(object=ms1$mod, data=ms1$dat, dataList = ms1$dat, 
            formula = Y ~ 1 + bsignal(x=X1, s=svals, df=3))


## sof/ scalar response with functional effect 
ms1a <- fitModel(formula =  Y ~ bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = Gaussian(), nygrid = 1)  ## scalar response as functional response with 1 observation 
testMethods(object=ms1a$mod, data=ms1a$dat, dataList = ms1a$dat, 
            formula = Y ~ bsignal(x=X1, s=svals, df=3))


## median sof/ scalar response with functional effect
ms2 <- fitModel(formula =  Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = Laplace(), nygrid = 1)   
testMethods(object=ms2$mod, data=ms2$dat, dataList = ms2$dat, 
            formula = Y ~ 1 + bsignal(x=X1, s=svals, df=3), family = Laplace())


###########################################################################################
###################### models with other family 

############# median regression: linear function-on-function,  \int x(s)beta(s,t)ds
m1b <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = Huber() )

m1b <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = ExpectReg() )

m1b <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = Laplace() )

m1b <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = QuantReg() )
## , numInt = "Riemann" does not work as family QuantReg() cannot deal with weights !!
# m1b$mod
testMethods(object=m1b$mod, data=m1b$dat, dataList = m1b$dat, 
            formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), family = QuantReg() )

m1b$mod$call 
## family in call is , family = ..1, does not happen if model is fitted directly in global environment
## because of this ugly family in call, the checks for update() do not work 
## TODO: find out why the family in call is incorrect 


############# gamma regression: linear function-on-function,  \int x(s)beta(s,t)ds
m1a <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = GammaReg() )
# m1a$mod
testMethods(object=m1a$mod, data=m1a$dat, dataList = m1a$dat, 
            formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), family = GammaReg())


############# binomial model: linear function-on-function,  \int x(s)beta(s,t)ds
m1c <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = Binomial(), offset = "scalar" )
testMethods(object=m1c$mod, data=m1c$dat, dataList = m1c$dat, 
            formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), family = Binomial())
## plotResiduals() is supposed not to work in this case

m1c$mod$call 

## TODO: validateFDboost() does not work with non-continuous response - repair that!

############# binomial model with AdaBoost: linear function-on-function,  \int x(s)beta(s,t)ds
m1d <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = AdaExp(), offset = "scalar" )
testMethods(object=m1d$mod, data=m1d$dat, dataList = m1d$dat, 
            formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), family = AdaExp())


############# Poisson model: linear function-on-function,  \int x(s)beta(s,t)ds
m1e <- fitModel(formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), scenario=c( "int", "ff"), 
                family = Poisson(), offset = "scalar" )
testMethods(object=m1e$mod, data=m1e$dat, dataList = m1e$dat, 
            formula=Y ~ 1 + bsignal(x=X1, s=svals, df=3), family = Poisson())


