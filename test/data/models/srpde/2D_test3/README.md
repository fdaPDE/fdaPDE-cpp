## 2D SRPDE - Test 3

| Type of model       | Non parametric SRPDE model with costant coefficients PDE regularization. |
|:--------------------|:-------------------------------------------------------------------------|
| Regularization      | Laplacian with costant diffusion tensor                                  |
| Covariates          | No                                                                       |
| Sampling            | Data sampled at mesh nodes                                               |
| Mesh                | unit square                                                              |
| Boundary conditions | No                                                                       |
| order               | 1                                                                        |
| smoothing parameter | 10.0                                                                     |

Possible R script generating the data

```r
x <- seq(0,1, length.out = 60)
y <- x
locations <- expand.grid(x,y)

mesh <- create.mesh.2D(locations)

# Test function
a1 <- 1
a2 <- 4
z <- function(p){
  a1*sin(2*pi*p[,1])*cos(2*pi*p[,2])+a2*sin(3*pi*p[,1])
}

# Exact solution (pointwise at nodes)
sol_exact <- z(locations)
ndati <- length(DatiEsatti)

# Add error to simulate data
set.seed(7893475)
ran <- range(DatiEsatti)
data <- DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))

PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)

options("scipen" = 999)
write.csv(format(data, digits=16), "z.csv")
```
