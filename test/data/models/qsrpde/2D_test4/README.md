## 2D SRPDE - Test 2

| Type of model       | Semiparametric SQRPDE model with laplacian regularization. |
|:--------------------|:-----------------------------------------------------------|
| Regularization      | Isotropic unitary laplacian                                |
| Covariates          | Yes                                                        |
| Sampling            | Areal                                                      |
| Mesh                | c shaped                                                   |
| Boundary conditions | No                                                         |
| order               | 1                                                          |
| smoothing parameter | 5.62341325190349110^(-3)                                   |

Possible R script generating the data

```r

## Load mesh
load("horseshoe2D_areal.RData")
nodes = mesh$nodes

mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)

# Read the correct incidence matrix
incidence_matrix = as.matrix(read.csv("incidence_matrix.csv")[, -1])
n = dim(incidence_matrix)[1]    # number of subdomains
nnodes = dim(nodes)[1]

# beta true
beta_true = 3

# Generate data and covariates

n = dim(incidence_matrix)[1]
integration_nodes = data.frame(matrix(nrow = n, ncol = 11))
names(integration_nodes) = c("T1","T2","N1","N2","N3","N4", "xmin", "xmax", "ymin", "ymax", "label")

tri = mesh$triangles
for(i in 1:n)
{
  tri_used = which(incidence_matrix[i,] == 1)
  integration_nodes$T1[i] = tri_used[1]
  integration_nodes$T2[i] = tri_used[2]
  nodes_used = unique(c(tri[tri_used[1],],tri[tri_used[2],]))
  integration_nodes$N1[i] = nodes_used[1]
  integration_nodes$N2[i] = nodes_used[2]
  integration_nodes$N3[i] = nodes_used[3]
  integration_nodes$N4[i] = nodes_used[4]
  integration_nodes$label[i] = i
  xvec = c(nodes[integration_nodes$N1[i],1],nodes[integration_nodes$N2[i],1],nodes[integration_nodes$N3[i],1],nodes[integration_nodes$N4[i],1])
  yvec = c(nodes[integration_nodes$N1[i],2],nodes[integration_nodes$N2[i],2],nodes[integration_nodes$N3[i],2],nodes[integration_nodes$N4[i],2])
  integration_nodes$xmin[i] = min(xvec)
  integration_nodes$xmax[i] = max(xvec)
  integration_nodes$ymin[i] = min(yvec)
  integration_nodes$ymax[i] = max(yvec)

}

area_domains = rep(0.04, n)

sol_integrand = matrix(nrow = n , ncol = 1)
for(i in 1:n){
  v = c(z(mesh$nodes[integration_nodes$N1[i],]), z(mesh$nodes[integration_nodes$N2[i],]), z(mesh$nodes[integration_nodes$N3[i],]), z(mesh$nodes[integration_nodes$N4[i],]))
  
  res = integral2(fun = z_xy, xmin = integration_nodes$xmin[i], xmax = integration_nodes$xmax[i] ,
                  ymin = integration_nodes$ymin[i], ymax = integration_nodes$ymax[i], 
                  vectorized = FALSE)
  sol_integrand[i] = (res$Q)/area_domains[i]
}

nnodes = dim(mesh$nodes)[1]
sol_exact <- numeric(nnodes)

set.seed(3)
covariates=matrix(0,nrow= n ,ncol=1)
covariates[,1]=rbeta(n ,shape1=2,shape2=2)  # sampling covariates from beta distr

write.csv(format(covariates, digits=16), "X.csv")

y_true <- sol_integrand + covariates%*%beta_true

# Add (skewed zero-mean) noise 
xi_ = 4 
omega_ = 0.1*(max(y_true)-min(y_true))
alpha_noise = 5
delta_ = alpha_noise / (sqrt(1+alpha_noise^2))
scale_noise = 1

noise_true_mean = scale_noise*(xi_ + omega_*delta_*sqrt(2/pi))
noise_true_sd = scale_noise*sqrt(omega_*(1-2*delta_^2/pi))

for(i in 1:nnodes){
  sol_exact[i] <- z(mesh$nodes[i,])
}

set.seed(441)
skewed_t <- scale_noise*rsn(n, xi = xi_, omega = omega_, alpha = alpha_noise) 

skewed_t <- skewed_t - noise_true_mean   # zero mean noise   
noise_true_quantile = scale_noise*qsn(alpha, xi = xi_, omega = omega_, alpha = alpha_noise) - noise_true_mean

# Generate simulated data 
data <- y_true + skewed_t
write.csv(format(data, digits=16), "z.csv")

```
