# TEST FOR MESH C 2D :

# Model: 
#   Family: Gamma
#   Covariate: Yes
#   Location: at mesh nodes 

#.rs.restartR()
rm(list=ls())
graphics.off()


# LOAD library ----------------------------------------

library(fdaPDE)

# FAMILY CHOICE ----------------------------------------
FAMILY = "gamma" 
# define the link and inverse link functions
link<-function(x){-1/x}
inv.link<-link
# BETA ------------------------------------------------

beta1=-2/5
beta2=3/10
betas_truth = c(beta1,beta2)
# lambda ---------------------------------------------


lambda = 10^seq(0,2,length.out = 20)
GCVFLAG=T
# scale param -----------------------------------------

scale.param = .1

# mesh reading:
data(horseshoe2D)

nodes = rbind(horseshoe2D$boundary_nodes, horseshoe2D$locations)
mesh = fdaPDE::create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh = mesh, minimum_angle = 30 )

FEMbasis = fdaPDE::create.FEM.basis(mesh)

# PLOT mesh ------
{
  plot(mesh, lwd=1.5, cex = 0.4)
}


# locations -----------------------------------------------
loc = mesh$nodes
nloc = dim(loc)[1]

# mesh with points
{
  plot(mesh, lwd=1, cex = 0.4)
  par(new=TRUE)
  plot(loc, pch = 19 , cex = 0.8, col = "#101875")
}

# True Field ---------------------
# 
fs.test=function (x, y, r0 = 0.1, r = 0.5, l = 3, b = 1, exclude = FALSE) 
{
  q <- pi * r/2
  a <- d <- x * 0
  ind <- x >= 0 & y > 0
  a[ind] <- q + x[ind]
  d[ind] <- y[ind] - r
  ind <- x >= 0 & y <= 0
  a[ind] <- -q - x[ind]
  d[ind] <- -r - y[ind]
  ind <- x < 0
  a[ind] <- -atan(y[ind]/x[ind]) * r
  d[ind] <- sqrt(x[ind]^2 + y[ind]^2) - r
  ind <- abs(d) > r - r0 | (x > l & (x - l)^2 + d^2 > (r - 
                                                         r0)^2)
  f <- a * b + d^2
  if (exclude) 
    f[ind] <- NA
  attr(f, "exclude") <- ind
  f
}
a = 8;
b=10 # 10

z <- function(p, a = 8, b = 10){
  
  -(1/a)*(fs.test(p[1],p[2]) +b)
  
}

# exact solution
sol_exact=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact[i] <- -(1/a)*(fs.test(loc[i,1],loc[i,2]) +b) ## truth
}

range(sol_exact)

# covariates ---------------------------------------------------------
set.seed(42)
desmat=matrix(0,nrow=length(loc[,1]),ncol=2)
desmat[,1]=rbeta(length(loc[,1]),shape1=1.5,shape2=2)  # sampling covariates from beta distr.
desmat[,2]=rbeta(length(loc[,1]),shape1=3,shape2=2)+1  # sampling covariates from beta distr.

# Plot covariates
{
  par(mfrow=c(1,2))
  plot(desmat[,1], response, xlab = "first covariate", ylab = "response")
  plot(desmat[,2], response, xlab = "second covariate", ylab = "response")
}

# Mean parameter  ---------------------------------------------------------
param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]

mu<-inv.link(param)

# Sampling gamma data --------------------------------------------
response <- rgamma(length(mesh$nodes[,1]), shape=mu/scale.param, scale=scale.param)


# Fitting data --------------------------------------------
output_CPP_exact <- smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat, GCV=GCVFLAG, GCVmethod = "Exact",
                                    lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL)
output_CPP_stoc <- smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat, GCV=GCVFLAG, GCVmethod = "Stochastic",
                                   lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL)



limsupvec = rbind(sol_exact+0.5, output_CPP$fn.eval[,output_CPP$bestlambda] + 0.5)
liminfvec = rbind(sol_exact-0.5, output_CPP$fn.eval[,output_CPP$bestlambda] - 0.5)
limval = c(max(limsupvec), min(liminfvec))

# mesh ------
{
  plot(mesh, lwd=1.5, cex = 0.4)
}

# mesh with points ------
{
  plot(mesh, lwd=1, cex = 0.4)
  par(new=TRUE)
  plot(loc, pch = 19 , cex = 0.8, col = "#101875")
}
# true field -----------
{
  nx<-250
  ny<-100 
  xvec <- seq(-1,4,length=nx)
  yvec<-seq(-1,1,length=ny)
  fsb <- list(fs.boundary())
  
  z_mat = matrix(ncol = length(yvec), nrow = length(xvec))
  for(i in 1:length(xvec)){
    for(j in 1:length(yvec)){
      z_mat[i,j] <- z(c(xvec[i],yvec[j]))
    }
  }
  image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
  lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
  contour(xvec,yvec,z_mat,add=TRUE)
}
# sampled data ----------
{
  library(RColorBrewer)
  pal = colorRampPalette(c("red","yellow"))
  pal = colorRampPalette(brewer.pal(6,"Blues"))
  plot_order = findInterval(response,sort(response))

  plot(mesh, lwd=1, cex = 0.5)
  par(new=TRUE)
  plot(loc[,1], loc[,2], pch = 19 , cex = 1, col = pal(nloc)[plot_order], xlab = "x", ylab = "y") # TRUE (1)
}
# mu  ------------------ 
{
nx<-250
ny<-100 
xvec <- seq(-1,4,length=nx)
yvec<-seq(-1,1,length=ny)
fsb <- list(fs.boundary())
covariate_evaluation = W_plot%*%betas_truth

z_mat = matrix(ncol = length(yvec), nrow = length(xvec))
for(i in 1:length(xvec)){
  for(j in 1:length(yvec)){
    z_mat[i,j] <- inv.link(z(c(xvec[i],yvec[j])))
  }
}

image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
contour(xvec,yvec,z_mat,add=TRUE)
}

# covariates  ---------
{
  par(mfrow=c(1,2))
  plot(desmat[,1], response, xlab = "first covariate", ylab = "response")
  plot(desmat[,2], response, xlab = "second covariate", ylab = "response")
}

# estimated field: exact GCV ----------

fem_exact <- FEM(output_CPP_exact$fit.FEM$coeff[,output_CPP_exact$bestlambda], FEMbasis)
#image(fem_exact)
nx<-250
ny<-100 
xvec <- seq(-1,4,length=nx)
yvec<-seq(-1,1,length=ny)
fsb <- list(fs.boundary())
z_mat <- matrix(ncol = length(yvec), nrow = length(xvec))
for(i in 1:length(xvec)){
  loc_vec = cbind(rep(xvec[i],length(yvec)),yvec)
  z_mat[i,] = eval.FEM(fem_exact, locations = loc_vec)
}

png("fitted_field_exactGCV_gamma2D(meshC).png", width = 720, height = 720, res = 150)
image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
contour(xvec,yvec,z_mat,add=TRUE)
dev.off()

# estimated field: Stochastic GCV ------

fem_stoc <- FEM(output_CPP_stoc$fit.FEM$coeff[,output_CPP_stoc$bestlambda], FEMbasis)
#image(fem_stoc)

z_mat <- matrix(ncol = length(yvec), nrow = length(xvec))
for(i in 1:length(xvec)){
  loc_vec = cbind(rep(xvec[i],length(yvec)),yvec)
  z_mat[i,] = eval.FEM(fem_stoc, locations = loc_vec)
}

image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
contour(xvec,yvec,z_mat,add=TRUE)


  
  