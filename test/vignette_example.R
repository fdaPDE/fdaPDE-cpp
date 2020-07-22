#.rs.restartR()
rm(list=ls())
graphics.off()
# LOAD library 
library(fdaPDE)


#### Smoothing with binomial response (GSR method) ####
data(binomialSquare2D)

FEMbasis = create.FEM.basis(square2D.mesh)

# 2D true field (function f)
a1=-2.5
a2=0.8

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1])
  
}

# exact solution
sol_exact=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact[i] <- z(loc[i,])
}
# solution on nodes
nnodes = dim(square2D.mesh$nodes)[1]
sol_nodes = numeric(nnodes)
for(i in 1:nnodes){
  sol_nodes[i] = z(square2D.mesh$nodes[i,])
}

# Choose lambda with GCV:
lambda = 10^seq(-5,0,length.out = 20)
GCVFLAG = TRUE
GCVmethod='Exact'

# Fitting solution for the binomial response
response # binomial samples (true or false)
FAMILY = "binomial"
solution = smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL, GCV=GCVFLAG, GCVmethod = GCVmethod,
                                 lambda = lambda, max.steps.FPIRLS = 15, family=FAMILY, mu0=NULL, scale.param=NULL)

plot(log10(lambda),solution$GCV)

best_lambda = solution$bestlambda
best_fit = solution$fit.FEM$coeff[,best_lambda]

image(FEM(sol_nodes, FEMbasis))
image(FEM(best_fit, FEMbasis))


