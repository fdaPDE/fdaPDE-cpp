# TEST for mesh HUB (2.5D)
# Model: 
#   -Family: "poisson"
#   - Covariates
#   - Locations at mesh nodes


#.rs.restartR()
rm(list=ls())
graphics.off()


# LOAD library ----------------------------------------
library(fdaPDE)
library(purrr)

# FAMILY CHOICE ----------------------------------------

FAMILY = "poisson"

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# BETA ------------------------------------------------

beta1= 2
beta2= -3
betas_truth = c(beta1,beta2)

# lambda ---------------------------------------------

lambda = c(0.0001,0.001,0.005,0.006,0.007,0.008,0.009,0.0095, 0.00975,0.01,0.0125,0.015,0.02,0.03,0.04,0.05,0.1,1,10,100,1000)
GCVFLAG=T

# scale param -----------------------------------------

scale.param = 1


# mesh reading: -----------------------------------------------
data(hub25Ddata)

mesh <- create.mesh.2.5D(nodes = nodes,triangles = triangles)

FEMbasis <- create.FEM.basis(mesh)

# locations -----------------------------------------------
loc=mesh$nodes
nloc = dim(loc)[1]

# 2.5D random field (function f) ------------------------------
a1 = rnorm(1,mean = 0, sd = 1)
a2 = rnorm(1,mean = 0, sd = 1)
a3 = rnorm(1,mean = 0, sd = 1)

sol_exact = numeric(nloc)
for (i in 0:(nloc-1)){
  sol_exact[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) + 7
}
range(sol_exact) 


# covariates ---------------------------------------------------------
set.seed(42)

desmat=matrix(0,nrow=nloc, ncol=2)

desmat[,1]= sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
desmat[,2]= rnorm(nloc, mean=2, sd=0.1)


# samoploing response: --------------------------------------------------  
ran=range(desmat%*%betas_truth + sol_exact) 
param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]

mu<-inv.link(param)
range(mu)
response <- rpois(nloc, lambda = mu)

# Fitting --------------------------------------------------------
output_CPP_exact<- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat, GCV=GCVFLAG, GCVmethod = "Exact",
                         lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL)
output_CPP_stoc <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat, GCV=GCVFLAG, GCVmethod = "Stochastic",
                                 lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL, scale.param=NULL)

beta_exact = output_CPP_exact$beta
beta_stoc = output_CPP_stoc$beta

func_estimation_exact = output_CPP_exact$fit.FEM$coeff[,output_CPP_exact$bestlambda]
func_estimation_stoc = output_CPP_stoc$fit.FEM$coeff[, output_CPP_stoc$bestlambda]


# Results: -----------------------------------------------------------------------------------------------------

plot(log10(lambda),output_CPP_exact$GCV)
plot(log10(lambda),output_CPP_stoc$GCV)

plot(FEM(sol_exact, FEMbasis)) # Plot the True Field
plot(FEM(func_estimation_exact,FEMbasis)) # Plot the estimated field using Exact GCV
plot(FEM(func_estimation_stoc,FEMbasis)) # Plot the estimated field using Stochastic GCV


