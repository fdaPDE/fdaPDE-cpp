library(fdaPDE2)

## load data
data(quasicircle2Dareal)
mesh = quasicircle2Dareal$mesh
incidence_matrix = quasicircle2Dareal$incidence_matrix

mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
)

## define regularizing term
pde <- new(SpaceVarying_2D_Order1, mesh_data)

## set null dirichlet boundary conditions
nodes <- pde$get_dofs_coordinates()
dirichletBC <- as.matrix(rep(0., times = dim(nodes)[1]))
pde$set_dirichlet_bc(dirichletBC)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## compute PDE parameters at quadrature nodes
n_quadrature_nodes <- length(quadrature_nodes)
## diffusion tensor
R <- 2.8; K1 <- 0.1; K2 <- 0.2
K_func <- function(p) { 
    10*rbind(
        c(p[2]^2 + K1 * p[1]^2 + K2 * (R^2 - p[1]^2 - p[2]^2), (K1 - 1) * p[1] * p[2]),
        c((K1 - 1) * p[1] * p[2], p[1]^2 + K1 * p[2]^2 + K2 * (R^2 - p[1]^2 - p[2]^2))
    )
}
## transport vector
beta <- 0.5
b_func <- function(p) { 
    10 * beta * c(p[1], p[2])
}

PDE_parameters <- list("diffusion" = as.matrix(t(apply(quadrature_nodes, 1, K_func))),
                       "transport" = as.matrix(t(apply(quadrature_nodes, 1, b_func))), 
                       "reaction"  = as.matrix(rep(0., times = n_quadrature_nodes)))
pde$set_PDE_parameters(PDE_parameters)

## define SRPDE model
model <- new(SRPDE_SpaceVarying_2D_Areal, pde)
model$set_subdomains(incidence_matrix)

## Add error to simulate data
set.seed(5839745)
sol_exact <- quasicircle2Dareal$data
data <- sol_exact + rnorm(length(sol_exact), sd = 0.05 * (max(sol_exact) - min(sol_exact)))

## set smoothing parameter
lambda <- 10^-3
model$set_lambda_s(lambda)
## set observations
model$set_observations(as.matrix(data))
## solve smoothing problem
model$solve()
