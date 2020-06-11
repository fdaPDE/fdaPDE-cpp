# TESTS FOR CHAROTID 2D

# Model: 
#   -Family: "poisson"
#   - NO Covariates
#   - Location: at mesh nodes ( in test 5.2.6 there are different locations)
#   - Space Varying ( forcing term u)
#   - Dirichlet Boudary conditions

# Simulation with the same True field and PDE of the test in section 5.2.6


#.rs.restartR()
rm(list=ls())
graphics.off()


# LOAD library ----------------------------------------

library(fdaPDE)
library(SuppDists)

# FAMILY CHOICE ----------------------------------------

FAMILY = "poisson"

  # load the poisson links function, in order to simulate the poisson data
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv

# lambda ---------------------------------------------

lambda = c(0.5, 0.75, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.25, 1.5, 2, 3, 4, 8, 10)
GCVFLAG=T

# scale param -----------------------------------------
scale.param = 1

# MESH ---------------------------

# mesh reading:
data(peak2Ddata)
# create and plot the triangulation
mesh <- create.mesh.2D(nodes, order = 1)
ri <- mesh$triangles
plot(mesh,asp=1)

FEMbasis <- fdaPDE::create.FEM.basis(mesh)

# locations -----------------------------------------------

loc = mesh$nodes

# 2D random field (function f) ------------------------------
z<-function(p)
{
  if(p[1] + p[2] == 0){
    0
  }else{
    2*(sin(sqrt(p[1]^2 + p[2]^2)*pi/2))
  }
}

# exact solution
sol_exact=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact[i] <- z(loc[i,])
}

range(sol_exact) 

# mean parameter ---------------------------------------------------------
param = sol_exact


range(param)

mu<-inv.link(param)

range(mu)

# PDE parameters -----------
{
  # anisotropy matrix
  sigma <- matrix(c(5/6,0,0,1/6),nrow=2,ncol=2,byrow=TRUE)
  
  K_func<-function(points)
  {
    output = array(0, c(2, 2, nrow(points)))
    for (i in 1:nrow(points))
      output[,,i] = sigma
    output
  }
  b_func<-function(points)
  {
    output = array(0, c(2, nrow(points)))
    for (i in 1:nrow(points))
      output[,i] = 0
    output
  }
  
  c_func<-function(points)
  {
    rep(c(0), nrow(points))
  }
  # forcing function
  u_func<-function(points)
  {
    output = array(0, c(1, nrow(points)))
    for (i in 1:nrow(points))
      output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,1,0)
    output
  }
  PDE_parameters <- list(K = K_func, b = b_func, c = c_func, u = u_func)
}

# Boundary Condition -----------------
{
  # homogeneous dirichlet boundary conditions
  BC <- list(BC_indices=which(nodes[,1]^2+nodes[,2]^2>24.5),
             BC_values=rep(0,length(which(nodes[,1]^2+nodes[,2]^2>24.5))))
}


# sampling poisson response: ---------------
set.seed(42)
response <- rpois(length(mu), lambda = mu)

  
# Fitting -------------

output_CPP_exact <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL, GCV=GCVFLAG, GCVmethod = "Exact", BC =BC, PDE_parameters = PDE_parameters, lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL)
  
output_CPP_stoc <- smooth.FEM(location = NULL, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL, GCV=GCVFLAG, GCVmethod = "Stochastic", BC =BC, PDE_parameters = PDE_parameters, lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL)

func_estimation_exact = output_CPP_exact$fit.FEM$coeff[, output_CPP_stoc$bestlambda]
func_estimation_stoc = output_CPP_stoc$fit.FEM$coeff[, output_CPP_stoc$bestlambda]
  
# Resutls: ----------
image(FEM(sol_exact, FEMbasis)) #true Field
image(FEM(func_estimation_exact,FEMbasis)) # estimation function with exact GCV
image(FEM(func_estimation_stoc,FEMbasis)) # #estimation function with stochastic GCV





