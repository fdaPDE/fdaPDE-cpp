## 2D SRPDE - Test 2

| Type of model       | Semiparametric SRPDE model with laplacian regularization. |
|:--------------------|:-----------------------------------------------------------|
| Regularization      | Isotropic unitary laplacian                                |
| Covariates          | Yes                                                        |
| Sampling            | Data sampled at given locations                            |
| Mesh                | c shaped                                                   |
| Boundary conditions | No                                                         |
| order               | 1                                                          |
| smoothing parameter | 0.2201047                                                  |

Possible R script generating the data

```r
data(horseshoe2D)
options("scipen" = 999)

mesh <- create.mesh.2D(
	   nodes    = horseshoe2D$boundary_nodes, 
	   segments = horseshoe2D$boundary_segments
	)
locations <- refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
write.csv(format(locations, digits=16), "locs.csv")

mesh <- refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

ndata <- nrow(locations)

# Create covariates
set.seed(509875)
cov1 <- rnorm(ndata, mean = 1, sd = 2)
cov2 <- sin(locations[,1])
write.csv(format(cbind(cov1,cov2), digits=16), "X.csv")

# Exact solution (pointwise at nodes)
sol_exact = fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2

# Add error to simulate data
set.seed(543663)
ran <- range(sol_exact)
data <- sol_exact + rnorm(length(sol_exact), mean=0, sd=0.05*abs(ran[2]-ran[1]))

write.csv(format(data, digits=16), "z.csv")
```
