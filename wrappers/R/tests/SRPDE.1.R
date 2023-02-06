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
pde <- new(Laplacian_2D_Order1, mesh_data)

## set null dirichlet boundary conditions
nodes <- pde$get_dofs_coordinates()
dirichletBC <- as.matrix(rep(0., times = dim(nodes)[1]))
pde$set_dirichlet_bc(dirichletBC)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define SRPDE model
model <- new(SRPDE_Laplacian_2D_GeoStatNodes, pde)

## test function
f <- function(x, y, z = 1) {
    coe <- function(x, y) 1 / 2 * sin(5 * pi * x) * exp(-x^2) + 1
    return(sin(2 * pi * (coe(y, 1) * x * cos(z - 2) - y * sin(z - 2))) *
           cos(2 * pi * (coe(y, 1) * x * cos(z - 2 + pi / 2) + coe(x, 1) *
                         y * sin((z - 2) * pi / 2))))
}

## exact solution (pointwise at nodes)
sol_exact <- f(mesh$nodes[, 1], mesh$nodes[, 2])
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
