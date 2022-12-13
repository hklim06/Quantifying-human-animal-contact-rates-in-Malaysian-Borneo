###############################################################################################################
########################### Example NHP spatial model using INLA  ############################
###############################################################################################################
rm(list=ls())

library(plyr)
library(sp)
library(INLA)
library(gridExtra)
library(lattice)
library(fields)
library(dplyr)
library(ggplot2)
library(ROCR)
library(rgdal)
library(rgeos)
library(raster)
library(MASS)
library(excursions)
library(devtools)
library(brinla)
library(raster)

## Data preparation

## Import your input file
df <- read.csv("/users/R/input_file.csv")
names(df)[2] <- "n"
names(df)[3] <- "y"
df$houseID <- NULL
df$link <- NA
df$X <- NULL
df$Sample <- NULL

## Load prediction data
pr <- read.csv("/users/R/prediction_layer.csv")
pr$X <- NULL
pr$y <- NA
pr$n <- 1
pr$link <- 1


## Merge data
df2 <- rbind(df, pr)


## Create unique ID for each place sampled
df2$xy <- paste(df2$GPS_XSS_X, df2$GPS_XSS_Y)
df2$sp <- match(df2$xy, unique(df2$xy))
names(df2)[2] <- "monkey"

## Scale environmental variables
en <- data.frame(scale(a[5:25]), center=TRUE, scale=TRUE)
en$center <- NULL; en$scale <- NULL


## Variable selection - stepAIC
X <- en
X$roads <- NULL
X$lc <- NULL
Y <- a[c("monkey", "n")]; colnames(Y) <- c("y", "n")
XY <- cbind(X, Y)
XY <- na.omit(XY)
XY$f <- XY$n - XY$y
XY$n <- NULL
mf <- glm(cbind(y,f)~ ., family="binomial", data=XY)
summary(mf)
ms <- mf %>% stepAIC(direction = "both")
ms


## Select variables for the model 
sp <- df2$sp #household ID
n <- df2$n #number of trials
Y <- df2$monkey #number of positive outcomes

inla.df <- data.frame(sp, en, n, Y)
covars <- inla.df[,c(1:22)] 


# Scale coordinates
xmin0 <- min(df2$GPS_XSS_X)
ymin0 <- min(df2$GPS_XSS_Y)
df2$x <- df2$GPS_XSS_X - xmin0
df2$y <- df2$GPS_XSS_Y - ymin0
df2$x <- df2$x/1000
df2$y <- df2$y/1000

## Set domain
xmin <- min(df2$x)
ymin <- min(df2$y)
xmax <- max(df2$x)
ymax <- max(df2$y)
x <- c(xmin, xmin, xmax, xmax)
y <- c(ymin, ymax, ymin, ymax)
boundary <- cbind(x,y)

## Create mesh
pts <- cbind(df2$x, df2$y)
mesh1 <- inla.mesh.create.helper(points=pts, max.edge=c(2,4), cut=1)
plot(mesh1)
points(pts[,1], pts[,2], pch=19, cex=.5, col="red")


## Create observation matrix
A <- inla.spde.make.A(mesh1, loc=pts)
ind <- inla.spde.make.index('s', mesh1$n)

# Priors: precision = 1/ sigma ^2 - weakly informative priors
fixed.priors <- list(mean.intercept = 0, prec.intercept=1/100, mean = list(ndvi=0, roads=0, pop=0, dem=0, slp=0, asp=0, sea=0, lc=0, avgt=0, mdr=0, maxt=0, prec=0, 
                                                                           seas=0, bf=0, mg=0, rb=0, ag=0, op=0, ir=0, of=0), 
                     prec=list(ndvi=1/10, roads=1/10, pop=1/10, dem=1/10, slp=1/10, asp=1/10, sea=1/10, lc=1/10, avgt=1/10, mdr=1/10, maxt=1/10, prec=1/10, 
                               seas=1/10, bf=1/10, mg=1/10, rb=1/10, ag=1/10, op=1/10, ir=1/10, of=1/10))

spde <- inla.spde2.matern(mesh1, alpha = 2)

## Create data stack
stk <- inla.stack(data=list(Y=Y), A=list(A,1), tag="sdata",
                  effects=list(ind, list(data.frame(b0=1, covars))))


############## MODEL FITTING - BINOMIAL ################

## Model 0: no spatial or random effect
form0 <- Y ~ 0 + b0 +ndvi + pop + dem + asp + sea + avgt + 
  mdr + maxt + mint + prec + seas + bf + mg + rb + ag + op + 
  ir + of
res0  <- inla(form0, family = "binomial", data=inla.stack.data(stk), Ntrials=n, 
              control.predictor = list(A=inla.stack.A(stk), compute=TRUE),
              control.fixed=fixed.priors, control.family=list(link="logit"),
              control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
              control.inla = list(strategy = "gaussian", int.strategy = "eb")) 
summary(res0)
res0$dic$dic 
save(res0, file = "nonspatialmodel.rda")

## Plot results
bri.fixed.plot(res0)

# Model fit
str(sdat <- inla.stack.index(stk, 'sdata')$data)
fitted.values0 <- res0$summary.fitted.values$mean[sdat]
fitted.values0 <- sapply(res0$marginals.linear.predictor[sdat], function(m)
  inla.emarginal(inla.link.invlogit, m))
observed.values <- Y/n
plot(fitted.values0, observed.values)
qplot(df$GPS_XSS_X, df$GPS_XSS_Y, colour=fitted.values0)+scale_color_gradient(low="blue", high="red")

## AUC
ROC_auc <- performance(prediction(fitted.values0,observed.values),"auc")
ROC_auc@y.values[[1]] # AUC




### Model 1: including spatial correlation

## Fit model with spatial correlation
form1 <- Y ~ 0 + b0 +ndvi + pop + dem + asp + sea + avgt + 
  mdr + maxt + mint + prec + seas + bf + mg + rb + ag + op + 
  ir + of+ f(s, model=spde)
res1 <- inla(form1, family = "binomial", data=inla.stack.data(stk), Ntrials=n,
             control.predictor = list(A=inla.stack.A(stk), compute=TRUE),
             control.fixed=fixed.priors, control.family=list(link="logit"),
             control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
             control.inla = list(strategy = "gaussian", int.strategy = "eb"))
summary(res1)
save(res1, file = "spatial_model.rda")
res1$dic$dic 
index.pred <- inla.stack.index(stk, "sdata")$data

df2$spde <- res1$summary.fitted.values[index.pred, "mean"]
df2$spde.sd <- res1$summary.fitted.values[index.pred, "sd"]
ggplot(df2, aes(x, y, color = spde)) + geom_point() + theme_bw() +  scale_color_viridis(option="rocket")
ggplot(df2, aes(x, y, color = spde.sd)) + geom_point() + theme_bw() +  scale_color_viridis(option="mako")



## Plot distributions
bri.hyperpar.plot(res1)
bri.fixed.plot(res1)

## Spatial range, range that it's uncorrelated
summary(res1$summary.hyperpar)
r.f <- inla.spde2.result(res1, "s", spde, do.transf = TRUE)
plot.default(r.f$marginals.range.nominal[[1]], type = "l", xlab = "Practical range", 
             ylab = "Density")
plot.default(r.f$marginals.variance.nominal[[1]], type = "l", xlab = expression(sigma[x]^2), 
             ylab = "Density")

# Calculate spatial range
spde.est <- inla.spde2.result(inla = res1, name = "s", spde = spde, do.transf = TRUE)
inla.zmarginal(spde.est$marginals.range.nominal[[1]])

# Assess model fit
str(sdat <- inla.stack.index(stk, 'sdata')$data)
fitted.values1 <- sapply(res1$marginals.linear.predictor[sdat], function(m)
  inla.emarginal(inla.link.invlogit, m))
plot(fitted.values1, observed.values)
cor.test(fitted.values1, observed.values, method = "pearson")
qplot(a$GPS_XSS_X, a$GPS_XSS_Y, colour=fitted.values1)+scale_color_gradient(low="blue", high="red")+ theme_bw()


# Assess model fit
str(sdat <- inla.stack.index(stk, 'sdata')$data)
fitted.values1 <- res1$summary.fitted.values$mean[sdat]
fitted.values1 <- sapply(res1$marginals.linear.predictor[sdat], function(m)
  inla.emarginal(inla.link.invlogit, m))

observed.values <- Y/n
plot(fitted.values1, observed.values)
cor.test(fitted.values1, observed.values, method = "pearson")
qplot(df2$GPS_XSS_X, df2$GPS_XSS_Y, colour=fitted.values1)+scale_color_gradient(low="blue", high="red")+ theme_bw()

library(pROC)
install.packages("ROCR")
library(ROCR)
## AUC
auc(observed.values, fitted.values1)





## Posterior samples, model testing, binomial probability
samp <- inla.posterior.sample(n=1000, res1)
pred.fun <- function(l, nsample=nrow(df2)){
  lp <- l$latent[1:nsample]
  lp <- inla.link.invlogit(lp) #between zero and one
  z <- rbinom(nsample, size=1, lp)
}
samples <- sapply(samp, pred.fun)
samples.mean <- apply(samples, 1, mean)
samples.sd <- apply(samples, 1, sd)

qplot(df2$x, df2$y, colour=samples.mean) + scale_color_viridis(option="turbo") + labs(color = "Probability") +
  theme_void()



## Estimate probability
pred.fun2 <- function(l, nsample=nrow(df2)){
  lp <- l$latent[1:nsample] ## extract the linear predictor in logit scale
  lp  <- inla.link.invlogit(lp)
}
samples <- sapply(samp, pred.fun2)
samples.mean <- apply(samples, 1, mean)
qplot(df2$x, df2$y, colour=samples.mean)

## Excursion probability - 10%
ep <- excursions.mc(samples, type=">", u=0.10, alpha=0.01)
qplot(df$x, df$y, colour=ep$rho)

## Write data for prediction
r <- data.frame(cbind(df$GPS_XSS_X, df$GPS_XSS_Y, samples.mean, ep$rho))
names(r) <- c("x", "y", "p", "exc")

# Save file
write.csv(r, "file_path.csv")

