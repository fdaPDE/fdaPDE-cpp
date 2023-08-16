## 2D SRPDE - Test 1

| Type of model       | Non parametric SRPDE model with laplacian regularization. |
|:--------------------|:----------------------------------------------------------|
| Regularization      | Isotropic unitary laplacian                               |
| Covariates          | No                                                        |
| Sampling            | Data sampled at mesh locations                            |
| Mesh                | unit square                                               |
| Boundary conditions | No                                                        |
| order               | 1                                                         |
| smoothing parameter | 5.623413e-05                                              |

Possible R script generating the data

```r
x <- seq(0,1, length.out = 60)
y <- x
locations <- expand.grid(x,y)

mesh <- create.mesh.2D(locations)

# Test function
f <- function(x, y, z <- 1){
  coe <- function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2) + 
  coe(x,1)*y*sin((z-2)*pi/2)))
}

# Exact solution (pointwise at nodes)
sol_exact <- f(mesh$nodes[,1], mesh$nodes[,2])

# Add error to simulate data
set.seed(7893475)
ran <- range(sol_exact)
data <- sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

options("scipen" = 999)
write.csv(format(data, digits=16), "z.csv")
```
