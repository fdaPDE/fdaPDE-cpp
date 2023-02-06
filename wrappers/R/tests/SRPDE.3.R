library(fdaPDE2)

## create mesh, unit square
x <- seq(0, 1, length.out = 20)
y <- x
locations <- expand.grid(x, y)
mesh <- fdaPDE::create.mesh.2D(locations)

mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
)

## define regularizing term
pde <- new(ConstantCoefficients_2D_Order1, mesh_data)

## set null dirichlet boundary conditions
nodes <- pde$get_dofs_coordinates()
dirichletBC <- as.matrix(rep(0., times = dim(nodes)[1]))
pde$set_dirichlet_bc(dirichletBC)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## set PDE parameters
PDE_parameters <- list("diffusion" = matrix(c(1, 0, 0, 4), nrow = 2),
                       "transport" = c(0, 0), 
                       "reaction"  = 0)
pde$set_PDE_parameters(PDE_parameters)

## define SRPDE model
model <- new(SRPDE_ConstantCoefficients_2D_GeoStatNodes, pde)

## test function
a1 <- 1; a2 <- 4
f <- function(p) {
    a1 * sin(2 * pi * p[, 1]) * cos(2 * pi * p[, 2]) + a2 * sin(3 * pi * p[, 1])
}
## exact solution (pointwise at nodes)
sol_exact <- f(locations)
## Add error to simulate data
set.seed(7893475)
ran  <- range(sol_exact)
data <- sol_exact + rnorm(dim(nodes)[1], mean = 0, sd = 0.05 * abs(ran[2] - ran[1]))

## set smoothing parameter
lambda <- 10^-3
model$set_lambda_s(lambda)
## set observations
model$set_observations(as.matrix(data))
## solve smoothing problem
model$solve()
