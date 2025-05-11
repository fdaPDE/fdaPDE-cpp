## 2D SRPDE - Test 4

| Type of model       | Non parametric SRPDE model with space-varying PDE regularization.                    |
|:--------------------|:-------------------------------------------------------------------------------------|
| Regularization      | non-constant coefficients diffusion-advection PDE with non homoegeneous forcing term |
| Covariates          | No                                                                                   |
| Sampling            | Areal                                                                                |
| Mesh                | quasi circle                                                                         |
| Boundary conditions | No                                                                                   |
| order               | 1                                                                                    |
| smoothing parameter | 0.001                                                                                |

Possible R script generating the data

```r
data(quasicircle2Dareal)
options("scipen" = 999)

mesh = quasicircle2Dareal$mesh
incidence_matrix = quasicircle2Dareal$incidence_matrix
write.csv(format(incidence_matrix, digits=16), "incidence_matrix.csv")

sol_exact = quasicircle2Dareal$data

# Add error to simulate data
set.seed(5839745)
data = sol_exact + rnorm(length(sol_exact), sd = 0.05*(max(sol_exact)-min(sol_exact)))
write.csv(format(data, digits=16), "z.csv")

# get quadrature nodes from the model
...

# Set space-varying PDE parameters
R = 2.8
K1 = 0.1
K2 = 0.2
beta = 0.5

K_func<-function(points){
    10*rbind(
        c(points[2]^2 + K1*points[1]^2 + K2*(R^2 - points[1]^2 - points[2]^2), 
          (K1-1)*points[1]*points[2]),
        c((K1-1)*points[1]*points[2], 
          points[1]^2 + K1*points[2]^2 + K2*(R^2 - points[1]^2 - points[2]^2))
    )
}

b_func<-function(points){
	10*beta*c(points[1],points[2])
}

u_func<-function(points){
	-ifelse((points[1]^2+points[i,2]^2) < 1, 100, 0)
}

diffusion_data = as.matrix(t(apply(quad_nodes, 1, K_func)))
write.csv(format(diffusion_data, digits=16), "K.csv")

advection_data = as.matrix(t(apply(quad_nodes, 1, b_func)))
write.csv(format(advection_data, digits=16), "b.csv")

forcing_data = as.matrix(apply(quad_nodes, 1, u_func))
write.csv(format(forcing_data, digits=16), "force.csv")
```
