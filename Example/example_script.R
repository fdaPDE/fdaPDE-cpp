# This script contains a few examples of new functionalities
# added in the context of the project.

# Change the following path to the current address 
# of the Example folder (i.e. the folder containing this script)
# From inside RStudio, one can also use:
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("~/Desktop/PACS/Example")

library(fdaPDE)

### New Mesh handling functionalities #####################
# 2.5D example
rm(list=ls())
load("./hub25Ddata.RData")
# Build mesh object (order 2 is also available)
mesh_hub<-create.mesh.2.5D(nodes=nodes,
                           triangles=triangles, 
                           order=1)
# Examine object contents
str(mesh_hub)
# mesh_hub is a list containing additional information (new):
# boundary markers (nodesmarkers and edgesmarkers), the edges matrix
# and the adjacency matrix (neighbors) are all new additions

# Plot mesh object
plot(mesh_hub)

# Extract boundary (new)
boundary_edges<-mesh_hub$edges[mesh_hub$edgesmarkers,]
boundary_nodes<-mesh_hub$nodes[mesh_hub$nodesmarkers,]
# Highlight boundary nodes
rgl.points(boundary_nodes[,1],
           boundary_nodes[,2],
           boundary_nodes[,3], 
           size=10, col="red")
# Highlight boundary edges
boundary_edges<-as.vector(t(boundary_edges))
rgl.lines(mesh_hub$nodes[boundary_edges,1],
          mesh_hub$nodes[boundary_edges,2],
          mesh_hub$nodes[boundary_edges,3], 
          lwd=5, col="red")

# Split mesh (new)
mesh_hub_split<-split(mesh_hub)
# Plot splitted mesh
plot(mesh_hub_split)

# 3D example
rm(list=ls())
load("./sphere3Ddata.RData")

# Build mesh object
# Order 2 is also available (new)
mesh_sphere<-create.mesh.3D(nodes=sphere3Ddata$nodes,
                            tetrahedrons=sphere3Ddata$tetrahedrons, 
                            order=1)
# Examine object contents
str(mesh_sphere)
# mesh_sphere is a list containing additional information (new):
# boundary markers (nodesmarkers and facesmarkers), the faces matrix
# and the adjacency matrix (neighbors) are all new additions

# Plot mesh object
plot(mesh_sphere)
# Plot mesh object with higher transparency to see interior nodes
plot(mesh_sphere, alpha=.25)

# Extract boundary (new)
boundary_faces<-mesh_sphere$faces[mesh_sphere$facesmarkers,]
boundary_nodes<-mesh_sphere$nodes[mesh_sphere$nodesmarkers,]
# Highlight boundary nodes
rgl.points(boundary_nodes[,1],boundary_nodes[,2],boundary_nodes[,3], size=10, col="red")
# Plot function already plots only boundary faces in this case

# Split mesh (new)
mesh_sphere_split<-split(mesh_sphere)
# Plot splitted mesh
plot(mesh_sphere_split)

### Second order FEM in 3D ################################
# Compare the accuracy and robustness of first 
# and second order methods in 3D settings.
# Order 2 in 3D is a new functionality.

rm(list=ls())

# Function to generate random points in a sphere
rsphere <- function(n, r = 1.0, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}

# Percentage of noise standard deviation on data range
# Try also different values, e.g 0.01 or 1
noisepercent<-0.25

# Vectors to contain RMSE evaluations
errors<-vector("list",2)
minerrors<-vector("numeric",2)
names(minerrors)<-c("Order 1", "Order 2")

for(j in 1:2){
  GCVFLAG=FALSE
  load("./sphere3Ddata.RData")
  
  # Build spherical mesh of order j (order 2 is new)
  mesh_sphere<-create.mesh.3D(sphere3Ddata$nodes,
                              sphere3Ddata$tetrahedrons,
                              order=j)
  
  FEMbasis <- create.FEM.basis(mesh_sphere)
  
  # Set smoothing parameters
  lambda= 10^seq(-6,1,by=.25)
  
  set.seed(5847947)
  
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  nnodes<-nrow(mesh_sphere$nodes)
  
  # Evaluate exact solution on mesh nodes
  func_evaluation = a1* sin(2*pi*mesh_sphere$nodes[,1]) +  a2* sin(2*pi*mesh_sphere$nodes[,2]) +  a3*sin(2*pi*mesh_sphere$nodes[,3]) + 1
  
  ran=range(func_evaluation)
  
  exact_sol=func_evaluation
  
  # Plot exact solution
  plot(FEM(exact_sol, FEMbasis))
  title3d(main=paste("Exact solution, order ", j), 
          col="black")
  
  # Add noise to exact solution to generate data for the simulation
  data= exact_sol + rnorm(nnodes,
                          mean=0,
                          sd=noisepercent * (ran[2]-ran[1]))
  
  # Compute the solution for each lambda
  output_CPP <- smooth.FEM(observations = data, 
                           FEMbasis = FEMbasis, 
                           lambda = lambda, 
                           GCV = GCVFLAG)
  
  # Generate locations inside the sphere
  nloc = 250
  loc=rsphere(nloc, 0.95)
  
  # Evaluate exact solution on generated locations
  func_evaluation_loc = a1* sin(2*pi*loc[,1]) +  a2* sin(2*pi*loc[,2]) +  a3*sin(2*pi*loc[,3]) + 1
  
  ran_loc=range(func_evaluation_loc)
  
  exact_sol_loc=func_evaluation_loc
  
  # Function to evaluate RMSE
  errorFun<-function(i){
    sqrt(sum((eval.FEM(output_CPP$fit.FEM, locations=loc, search="tree")[,i] - exact_sol_loc)^2)/nloc)
  }
  
  # Compute RMSE for each lambda
  errors[[j]]<-vapply(1:length(lambda),
                      errorFun,
                      FUN.VALUE = numeric(1))
  
  # Compute minimum RMSE (equivalently, find optimal lambda)
  minerrors[j]<-min(errors[[j]])
  
  # Plot optimal solution
  plot(FEM(coeff=output_CPP$fit.FEM$coeff[,which.min(errors[[j]])], 
           FEMbasis = FEMbasis))
  title3d(main=paste("Result, optimal lambda, order ", j), 
          col="black")
  
}

# Plot RMSE with respect to lambda
plot(log10(lambda), errors[[1]], 
     ylim=range(errors[[1]],errors[[2]]), 
     ylab="RMSE",
     main="3D 2nd Order Accuracy comparison")
points(log10(lambda), errors[[2]], 
       col="red")
legend("bottomright", legend=c("Order 1", "Order 2"),
       fill=c("black", "red"),cex=0.8)

# Print optimal RMSE
minerrors

### PDE penalization in 3D ################################
# Include PDE parameters and use a penalization term with
# anisotropic diffusion and possibly advection and reaction
# terms. All this is a new functionality.

rm(list=ls())

# Function to generate random points in a sphere
rsphere <- function(n, r = 1.0, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}

GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic' #for stochastic GCV (default)

# Build mesh: Sphere
load("./sphere3Ddata.RData")

mesh_sphere<-create.mesh.3D(sphere3Ddata$nodes,
                            sphere3Ddata$tetrahedrons,
                            order=1)

FEMbasis <- create.FEM.basis(mesh_sphere)

set.seed(5847947)

# Exact test function
nnodes = nrow(mesh_sphere$nodes)

# Set smoothing parameter
lambda=10^seq(-6, 1, by=.25)

# Set PDE parameters (in this case they are constant)
PDE_parameters_anys = list(K = diag(c(1,.5,1)), b = c(0,0,0), c = -4*pi^2)

# Evaluate exact solution on mesh nodes
exact_sol =  sin(2*pi*mesh_sphere$nodes[,1]) +  2 * sin(2*pi*mesh_sphere$nodes[,2]) +  sin(2*pi*mesh_sphere$nodes[,3])

# Plot exact solution
plot(FEM(exact_sol,FEMbasis))
title3d(main="Exact solution", 
        col="black")

# Add noise to generate data - 10% level of noise
data=exact_sol + rnorm(nrow(mesh_sphere$nodes), mean=0, sd=0.10*diff(range(exact_sol)))

# Compute the solution for each lambda
output_CPP <- smooth.FEM(observations = data, PDE_parameters = PDE_parameters_anys,
                                                 FEMbasis = FEMbasis, lambda = lambda, GCV = GCVFLAG)
  
# Generate locations in the sphere
nloc = 250
loc=rsphere(nloc, 0.95)

# Evaluate true solution on generated locations
exact_sol_loc =  sin(2*pi*loc[,1]) +  2 * sin(2*pi*loc[,2]) +  sin(2*pi*loc[,3])
  
# Function to evaluate RMSE - here we use walking search. 
# Walking search is also a new addition for 3D data
errorFun<-function(i){
  sqrt(sum((eval.FEM(output_CPP$fit.FEM, locations=loc, search="walking")[,i] - exact_sol_loc)^2)/nloc)
}

# Compute RMSE for each lambda
errors<-vapply(1:length(lambda),
                    errorFun,
                    FUN.VALUE = numeric(1))

# Plot RMSE with respect to lambda
plot(log10(lambda), errors, 
     ylab="RMSE", 
     main="3D PDE Regression example")

# Plot optimal solution
plot(FEM(coeff=output_CPP$fit.FEM$coeff[,which.min(errors)], 
         FEMbasis = FEMbasis))
title3d(main="Result, optimal lambda", 
        col="black")

