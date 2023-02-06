library(fdaPDE2)

## load mesh data and create mesh object
data(horseshoe2D)
mesh <- fdaPDE::create.mesh.2D(
                    nodes    = horseshoe2D$boundary_nodes,
                    segments = horseshoe2D$boundary_segments
                )
locations <- fdaPDE::refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh <- fdaPDE::refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

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
model <- new(SRPDE_Laplacian_2D_GeoStatLocations, pde)
model$set_locations(locations)

## generate data
ndata <- nrow(locations)

# Create covariates
set.seed(509875)
cov1 <- rnorm(ndata, mean = 1, sd = 2)
cov2 <- sin(locations[, 1])

# Exact solution (pointwise at nodes)
sol_exact <- fdaPDE::fs.test(locations[, 1], locations[, 2]) + 2 * cov1 - cov2
# Add error to simulate data
set.seed(543663)
ran <- range(sol_exact)
data <- sol_exact + rnorm(length(sol_exact), mean = 0, sd = 0.05 * abs(ran[2] - ran[1]))

## set smoothing parameter
lambda <- 10^-3
model$set_lambda_s(lambda)
## set observations and covariates
model$set_observations(as.matrix(data))
model$set_covariates(cbind(cov1, cov2))
## solve smoothing problem
model$solve()
