## 2D SQRPDE - Test 3

| Type of model       | Non parametric SRPDE model with costant coefficients PDE regularization. |
|:--------------------|:-------------------------------------------------------------------------|
| Regularization      | Laplacian with costant diffusion tensor                                  |
| Covariates          | No                                                                       |
| Sampling            | Data sampled at mesh nodes                                               |
| Mesh                | unit square                                                              |
| Boundary conditions | No                                                                       |
| order               | 1                                                                        |
| smoothing parameter | 5.623413251903491*e-04                                                   |

Possible R script generating the data

```r
nxx = 21
x <- seq(0, 1, length.out = nxx)
y <- x
mesh <- fdaPDE::create.mesh.2D(expand.grid(x, y))
nodes = mesh$nodes
nnodes = dim(nodes)[1]

mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)

pde_data <- list(
  "diffusion" = cbind(c(1,0), c(0,4)), 
  "transport" = rep(0,2), 
  "reaction" = 0 
)

locations = nodes

## data

# utility for data generation 
data_creation_skewed = function(data_generation, alpha, locations, nodes, 
                                beta_true = NULL, covariates = NULL, seed){
  
  # skewed data generated as a deterministic function of the position + skewed-gaussian noise 
  #         y = data_generation(x,y) + rsn 
  
  n = dim(locations)[1]
  y_true <- data_generation(locations[,1], locations[,2]) 
  if(!is.null(covariates)) {
    q = dim(covariates)[2]
    if(is.null(beta_true)){
      beta_true = seq(1,q)
    }
    y_true = y_true + covariates%*%beta_true  
  }
  
  # Add (skewed zero-mean) noise 
  library(sn)  # Skewed normal Distribution
  xi_ = 4 
  omega_ = 0.05*(max(y_true)-min(y_true))
  alpha_noise = 5
  delta_ = alpha_noise / (sqrt(1+alpha_noise^2))
  scale_noise = 1
  
  set.seed(seed)
  skewed_t <- scale_noise*rsn(n, xi = xi_, omega = omega_, alpha = alpha_noise) 
  
  noise_true_mean = scale_noise*(xi_ + omega_*delta_*sqrt(2/pi))
  noise_true_sd = scale_noise*sqrt(omega_*(1-2*delta_^2/pi))
  
  skewed_t <- skewed_t - noise_true_mean   # zero mean noise   
  noise_true_quantile = scale_noise*qsn(alpha, xi = xi_, omega = omega_, alpha = alpha_noise) - noise_true_mean
  
  # Generate simulated data (pointwise at nodes) 
  data <- y_true + skewed_t
  
  # True field
  fn_true <- noise_true_quantile + data_generation(locations[,1], locations[,2]) 
  f_true <- noise_true_quantile + data_generation(nodes[,1], nodes[,2]) 
  
  return(list(data = data, beta_true = beta_true, fn_true = fn_true, f_true = f_true, 
              noise_true_mean = 0, noise_true_sd = noise_true_sd))
}  


data_generation <- function(x, y, z = 1){
  a1 <- 1
  a2 <- 4
  a1*sin(2*pi*(x+0.17))*cos(2*pi*y)+a2*sin(3*pi*(x+0.17))
}   
data_result = data_creation_skewed(data_generation = data_generation,
                                   alpha = alpha, locations = locations, nodes = nodes, 
                                   seed = seed)
data = data_result$data
fn_true = data_result$fn_true       # true field at locations 
f_true = data_result$f_true         # true field at mesh nodes 


options("scipen" = 999)
write.csv(format(data, digits=16), "z.csv")

```
