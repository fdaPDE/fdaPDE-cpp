# TEST SPHERE 3D
# Model: 
#   -Family: "gamma"
#   - Covariates
#   - Locations not at mesh nodes


#.rs.restartR()
rm(list=ls())
graphics.off()

# LOAD library ----------------------------------------

library(fdaPDE)
library(purrr)
#library(fda)

# FAMILY CHOICE ----------------------------------------

FAMILY = "gamma"

link<-function(x){-1/x}
inv.link<-link

# BETA ------------------------------------------------

beta1=0.8
beta2=1.2
beta_exact=c(beta1,beta2)

# lambda ---------------------------------------------

lambda = 10^seq(-1,3,length.out = 10)
GCVFLAG=T
GCVmethod = 'Stochastic'
GCVmethod = 'Exact'

# scale param -----------------------------------------

scale.param = 1

# mesh reading:
data("sphere3Ddata")
sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
#plot(sphere3D)
FEMbasis <- fdaPDE::create.FEM.basis(sphere3D)
nodesLocations=sphere3D$nodes


set.seed(5847947)
set.seed(42)


# add evaluation points

new_grid3d = cbind(runif(1000,min = -1 , max = 1),runif(1000,min = -1 , max = 1),runif(1000,min = -1 , max = 1))
for(i in 1:length(new_grid3d[,1])){
  if(new_grid3d[i,1]^2 + new_grid3d[i,2]^2 + new_grid3d[i,3]^2 >= 0.98){
    new_grid3d[i,] <- c(NA,NA,NA)
  }
}

new_loc <- matrix(ncol = 3, nrow = sum(!is.na(new_grid3d[,1])))
new_loc[,1] <- new_grid3d[!is.na(new_grid3d[,1]),1] 
new_loc[,2] <- new_grid3d[!is.na(new_grid3d[,1]),2] 
new_loc[,3] <- new_grid3d[!is.na(new_grid3d[,1]),3] 

loc <- rbind(nodesLocations, new_loc)
dim(loc)

# Exact test function - locations at nodes
nnodes = sphere3D$nnodes
nloc = dim(loc)[1]
a1 = rnorm(1,mean = -2, sd = 1/8)
a2 = rnorm(1,mean = -2, sd = 1/8)
a3 = rnorm(1,mean = -2, sd = 1/8)

a1 = -1
a2 = -2
a3 = -3
func_evaluation = numeric(nloc)

for (i in 0:(nloc-1)){
  func_evaluation[i+1] = a1* sin(loc[i+1,1]) +  a2* sin(loc[i+1,2]) +  a3*sin(loc[i+1,3]) - 7
}
ran=range(func_evaluation) 
ran

# covariates

cov1_nonod=sin(2*pi*loc[,1])+sin((2*pi*loc[,2])^2)
cov2_nonod=cos(-2*pi*loc[,3])

plot(FEM(cov1_nonod[1:nnodes], FEMbasis))
plot(FEM(cov2_nonod[1:nnodes], FEMbasis))

W2=cbind(cov1_nonod,cov2_nonod)


# Evaluation on a planar cut
{
  grid_planar_ = cbind(runif(1000,min = -1 , max = 1),runif(1000,min = -1 , max = 1))
  radius = 1
  for(i in 1:length(grid_planar_[,1])){
    
    if(grid_planar_[i,1]^2+grid_planar_[i,2]^2 >= radius){
      grid_planar_[i,] <- c(NA,NA)
    }
    
  }    
  
  grid_planar <- matrix(ncol=2, nrow = sum(!is.na(grid_planar_[,1])))
  grid_planar[,1] <- grid_planar_[!is.na(grid_planar_[,1]),1]
  grid_planar[,2] <- grid_planar_[!is.na(grid_planar_[,2]),2]
  
  
  cov1_planar=sin(2*pi*grid_planar[,1])+sin((2*pi*grid_planar[,2])^2)
  meshgrid2d <- create.mesh.2D(nodes=grid_planar)
  plot(meshgrid2d)
  FEM_basis_planar = fdaPDE::create.FEM.basis(meshgrid2d)
  
  image(FEM(cov1_planar, FEM_basis_planar))
}


range(W2%*%beta_exact + func_evaluation) 

 

# data generation

theta = func_evaluation + W2%*%beta_exact
plot(FEM(theta[1:nnodes], FEMbasis))

mu = inv.link(theta)
plot(FEM(mu[1:nnodes], FEMbasis))


range(mu)

response <- rgamma(length(loc[,1]), shape=mu/scale.param, scale=scale.param)

plot(FEM(response[1:nnodes],FEMbasis))


output_CPP <- fdaPDE::smooth.FEM(locations = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = W2, GCV=GCVFLAG, GCVmethod = GCVmethod,
                                  lambda = lambda, max.steps.FPIRLS=15, family=FAMILY, mu0=NULL, scale.param=NULL)

plot(log10(lambda),output_CPP$GCV)

func_estimation = output_CPP$fit.FEM$coeff[,output_CPP$bestlambda]

plot(FEM(func_evaluation[1:nnodes],FEMbasis))
plot(FEM(func_estimation,FEMbasis))



