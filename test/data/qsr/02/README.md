## 2D SRPDE - Test 2

| Type of model       | Semiparametric SQRPDE model with laplacian regularization. |
|:--------------------|:-----------------------------------------------------------|
| Regularization      | Isotropic unitary laplacian                                |
| Covariates          | Yes                                                        |
| Sampling            | Data sampled at given locations                            |
| Mesh                | c shaped                                                   |
| Boundary conditions | No                                                         |
| order               | 1                                                          |
| smoothing parameter | 3.162277660168379*10^(-4)                                  |

Possible R script generating the data

```r
data(horseshoe2D)
options("scipen" = 999)

mesh_coarse <- create.mesh.2D(
  nodes    = horseshoe2D$boundary_nodes, 
  segments = horseshoe2D$boundary_segments
)

mesh <- refine.mesh.2D(mesh_coarse, maximum_area = 0.025, minimum_angle = 30)

## generate locations
locations <- refine.mesh.2D(mesh_coarse, maximum_area = 0.04)$nodes
write.csv(format(locations, digits=16), "locs.csv")

n <- nrow(locations)


## utility for data generation 
mean.function = function(locations){ 
  f = fs.test(locations[,1], locations[,2], b=1, exclude = FALSE)
  return(f)
} 

std.function = function(locations){ 
  f = fs.test(locations[,1], locations[,2], b=0.5, exclude = FALSE)
  min = min(f)
  f = 0.2*(f - 1.2*min(f))
  return(f)
}


## beta true 
beta_true = c(2, -1)

# Covariates
set.seed(21)
cov1 <- (rnorm(n, mean = 0 , sd = 1))
set.seed(55)
cov2 <- (rexp(n, rate = 1))
X = cbind(cov1 ,cov2)
write.csv(format(X, digits=16), "X.csv")

# Data
set.seed(543663)
data = rnorm(n, mean = mean.function(locations), 
             sd = std.function(locations)) + X%*%beta_true 

write.csv(format(data, digits=16), "z.csv")

```
