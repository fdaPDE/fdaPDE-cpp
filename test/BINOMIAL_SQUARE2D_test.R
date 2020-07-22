# BINOMIAL 2D SQUARE SIMULATION
# see section 5.2.1 on report for the theoretical setting

rm(list=ls())
graphics.off()

# LOAD library ----------------------------------------

library(fdaPDE)
library(purrr)

# LOAD data ----------------------------------------
data("BINOMIAL_SQUARE_2D")


# Description of binomial data variables -------

  mesh # 2D Square mesh
  loc # points over the mesh
  nloc # total number of points over the mesh
  mu # mean parameter
  response # binomial samples (true or false)
  sol_exact # true field (function z with coefficient a1 a2) valutated over the locations
  param # \theta(beta, sol_exact) in this example beta = NULL

# Plotting the Simulation Settings -----
# mesh 
{
  plot(mesh, lwd=1.5, cex = 0.4)
}
  
# mesh with points 
{
  plot(mesh, lwd=1, cex = 0.4)
  par(new=TRUE)
  plot(loc, pch = 19 , cex = 0.8, col = "#101875")
}
  
# true field
{
  xvec = seq(0,1,by = 0.01)
  yvec = seq(0,1,by = 0.01)
  xy = cbind(xvec,yvec)
  s = sort(xvec, index.return = TRUE)
  xvec = unique(xvec[s$ix])
  s = sort(yvec, index.return = TRUE)
  yvec = unique(yvec[s$ix])
  
  z_mat = matrix(ncol = length(yvec), nrow = length(xvec))
  for(i in 1:length(xvec)){
    for(j in 1:length(yvec)){
      z_mat[i,j] <- z(c(xvec[i],yvec[j]))
    }
  }

  image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
  contour(xvec,yvec,z_mat,add=TRUE)
}
  
# sampled data 
{
  plot(mesh, lwd=1, cex = 0.5)
  par(new=TRUE)
  plot(loc[response,], pch = 1 , cex = 1, col = "#0AA398") # TRUE (1)
  par(new=TRUE)
  plot(loc[!response,], pch = 4, cex = 1, col = "#E5401D" , xaxt = 'n', yaxt = 'n', ann=FALSE) # FALSE (0)
}
  
# mu (no covariates) 
{  
  xvec = seq(0,1,by = 0.01)
  yvec = seq(0,1,by = 0.01)
  xy = cbind(xvec,yvec)
  s = sort(xvec, index.return = TRUE)
  xvec = unique(xvec[s$ix])
  s = sort(yvec, index.return = TRUE)
  yvec = unique(yvec[s$ix])
  
  z_mat = matrix(ncol = length(yvec), nrow = length(xvec))
  for(i in 1:length(xvec)){
    for(j in 1:length(yvec)){
      z_mat[i,j] <- inv.link(z(c(xvec[i],yvec[j])))
    }
  }
  
  image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
  contour(xvec,yvec,z_mat,add=TRUE)
}

# Fitting solution --------
FAMILY = "binomial"  
lambda = 10^seq(-5,0,length.out = 20)  

output_CPP_exact <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL, GCV=T, GCVmethod = "Exact",
                                       lambda = lambda, max.steps.FPIRLS=15, family=FAMILY, mu0=NULL, scale.param=NULL)

output_CPP_stoch <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL, GCV=T, GCVmethod = "Stochastic",
                                       lambda = lambda, max.steps.FPIRLS=15, family=FAMILY, mu0=NULL, scale.param=NULL)


# plot the best lambda which minimize the GCV 
plot(log10(lambda),output_CPP_exact$GCV)

# save the best lambda 
best_lambda_ex = output_CPP_exact$bestlambda
best_lambda_stoch = output_CPP_stoch$bestlambda

# save the best estimated function for both exact and stochastic GCV method
f_est_ex = output_CPP_exact$fit.FEM$coeff[,best_lambda_ex]
f_est_stoc = output_CPP_exact$fit.FEM$coeff[,best_lambda_stoch]

# ------------- Results plot --------------

# plot settings
{
  limsupvec = rbind(sol_exact+0.5, output_CPP$fn.eval[,output_CPP$bestlambda] + 0.5)
  liminfvec = rbind(sol_exact-0.5, output_CPP$fn.eval[,output_CPP$bestlambda] - 0.5)
  limval = c(max(limsupvec), min(liminfvec))
}
# estimated field: exact GCV ----------
{
  fem_exact <- FEM(f_est_ex, FEMbasis)
  
  
  z_mat <- matrix(ncol = length(yvec), nrow = length(xvec))
  for(i in 1:length(xvec)){
    loc_vec = cbind(rep(xvec[i],length(yvec)),yvec)
    z_mat[i,] = eval.FEM(fem_exact, locations = loc_vec)
  }
}
image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
contour(xvec,yvec,z_mat,add=TRUE)
# estimated field: Stochastic GCV ------
{
  fem_stoc <- FEM(f_est_stoc, FEMbasis)
  #image(fem_stoc)
  
  z_mat <- matrix(ncol = length(yvec), nrow = length(xvec))
  for(i in 1:length(xvec)){
    loc_vec = cbind(rep(xvec[i],length(yvec)),yvec)
    z_mat[i,] = eval.FEM(fem_stoc, locations = loc_vec)
  }

}
image(xvec,yvec,z_mat,col=heat.colors(100),xlab="x",ylab="y")
contour(xvec,yvec,z_mat,add=TRUE)


