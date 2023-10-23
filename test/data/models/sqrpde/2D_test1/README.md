## 2D SQRPDE - Test 1

| Type of model       | Non parametric SQRPDE model with laplacian regularization.|
|:--------------------|:----------------------------------------------------------|
| Regularization      | Isotropic unitary laplacian                               |
| Covariates          | No                                                        |
| Sampling            | Data sampled at mesh locations                            |
| Mesh                | unit square                                               |
| Boundary conditions | No                                                        |
| order               | 1                                                         |
| smoothing parameter | 1.778279*e-05                                             |

Possible R script generating the data

```r
nxx = 21 
x <- seq(0, 1, length.out = nxx)
y <- x
locations <- expand.grid(x, y)
n = dim(locations)[1]
mesh <- fdaPDE::create.mesh.2D(locations)
nodes = mesh$nodes

## data

# utility for data generation 
data_creation_hetero_grf = function(nxx, nodes, locations, covariates = NULL, beta_true = NULL, seed, 
                                    tau2_mean = 1, nu_mean = 2, rho_mean = 0.3, tau2_sd = 0.3, nu_sd = 2, 
                                    rho_sd = 0.6, anisotropy = NULL){
  
  # heteroscedastic data starting from true mean & std fields Gaussian Random Field
  #         y = rnorm(n, mean = field_mean, sd = field_std)
  
  library(geoR)
  
  seed.fixed = 10 # a fixed seed for the grf (we do not want that mu,sigma change with 
  # the simulation seed)
  
  n = dim(locations)[1]
  nnodes = dim(nodes)[1]
  
  if(!is.null(anisotropy)){    # anisotropic case
    aniso.pars = c(anisotropy$angle, anisotropy$intensity)
  } else{                      # isotropic case
    aniso.pars = c(0, 1) 
  }
  
  # Mean field
  rho_mean = rho_mean / (2*sqrt(nu_mean))
  cov.pars_mean = c(tau2_mean, rho_mean)
  
  set.seed(seed.fixed)
  mean_loc = grf(n, grid = locations, sqrt(n), sqrt(n), xlims = c(0, 1), ylims = c(0, 1),
                 nsim = 1, cov.model = "matern",
                 cov.pars = cov.pars_mean , 
                 kappa = nu_mean, 
                 aniso.pars = aniso.pars,
                 messages = FALSE)$data   # nugget = 0, lambda = 1, aniso.pars,
  set.seed(seed.fixed)
  mean_nodes = grf(nnodes, grid = nodes, sqrt(nnodes), sqrt(nnodes), xlims = c(0, 1), ylims = c(0, 1),
                   nsim = 1, cov.model = "matern",
                   cov.pars = cov.pars_mean , 
                   kappa = nu_mean, 
                   aniso.pars = aniso.pars,
                   messages = FALSE)$data   # nugget = 0, lambda = 1, aniso.pars,
  
  
  # Std field
  rho_sd = rho_sd / (2*sqrt(nu_sd))
  cov.pars_sd = c(tau2_sd, rho_sd)
  
  set.seed(seed.fixed)
  sd_loc = grf(n, grid = locations, sqrt(n), sqrt(n), xlims = c(0, 1), ylims = c(0, 1),
               nsim = 1, cov.model = "matern",
               cov.pars = cov.pars_sd , 
               kappa = nu_sd, 
               aniso.pars = aniso.pars, 
               messages = FALSE)$data  # nugget = 0, lambda = 1, aniso.pars,
  set.seed(seed.fixed)
  sd_nodes = grf(nnodes, grid = nodes, sqrt(nnodes), sqrt(nnodes), xlims = c(0, 1), ylims = c(0, 1),
                 nsim = 1, cov.model = "matern",
                 cov.pars = cov.pars_sd , 
                 kappa = nu_sd, 
                 aniso.pars = aniso.pars, 
                 messages = FALSE)$data  # nugget = 0, lambda = 1, aniso.pars,
  
  if(min(sd_loc) <= 0){
    cat("\n Adjusting negative values of the std gaussian field... \n")
    sd_loc[which(sd_loc <= 0)] = sd_loc[which(sd_loc <= 0)] - 1.2*sd_loc[which(sd_loc <= 0)]
  }
  if(min(sd_nodes) <= 0){  # if there are negative values 
    cat("\n Adjusting negative values of the std gaussian field... \n")
    sd_nodes[which(sd_nodes <= 0)] = sd_nodes[which(sd_nodes <= 0)] - 1.2*sd_nodes[which(sd_nodes <= 0)]
  }
  
  set.seed(seed)   # this seed is simulation dependent 
  data = rnorm(n, mean = mean_loc, sd = sd_loc)
  if(!is.null(covariates)){
    q = dim(covariates)[2]
    if(is.null(beta_true)){
      beta_true = seq(1,q)
    }
    data = data + covariates%*%beta_true  
  }
  
  # Field true
  fn_true <- mean_loc + qnorm(alpha) * sd_loc
  f_true <- mean_nodes + qnorm(alpha) * sd_nodes
  
  return(list(data = data, beta_true = beta_true, 
              mean.field_nodes = mean_nodes, std.field_nodes = sd_nodes,
              mean.field_loc = mean_loc, std.field_loc = sd_loc, 
              fn_true = fn_true, f_true = f_true))
} 

data_result = data_creation_hetero_grf(nxx = nxx, nodes = nodes, 
                                       locations = locations, seed = seed)
data = data_result$data
fn_true = data_result$fn_true       # true field at locations 
f_true = data_result$f_true         # true field at mesh nodes 
mean_true = data_result$mean.field  # true mean field
sd_true = data_result$std.field     # true std deviation field 

## save data
options("scipen" = 999)
write.csv(format(data, digits=16), "z.csv")

```
