# TESTING SQUARE 2D - exponential

#.rs.restartR()
rm(list=ls())
graphics.off()

# LOAD library ----------------------------------------

library(fdaPDE)

# FAMILY CHOICE ----------------------------------------

FAMILY = "exponential"  # "cloglog", "probit", "poisson", "binomial", "gamma", "exponential" , ( "gaussian", "inv_gaussian" )

link<-function(x){-1/x}
inv.link<-link 


# lambda ---------------------------------------------

lambda = 10^seq(-5,0,length.out = 20)
GCVFLAG=T
GCVmethod = 'Stochastic'
GCVmethod='Exact'

# scale param -----------------------------------------

scale.param = 1

# -------------------- CODE ---------------------------

# mesh reading:
data(square2Ddata)

mesh = fdaPDE::create.mesh.2D(nodes=nodes)
# x11()
plot(mesh, lwd=3, cex = 1.9)
# axis(1)
# axis(2)
nnodes = dim(mesh$nodes)[1]
FEMbasis = fdaPDE::create.FEM.basis(mesh)


# locations -----------------------------------------------
set.seed(42)
# locations
nloc = 800
xobs=runif(min=0,max=1,n=nloc)
yobs=runif(min=0,max=1,n=nloc)
loc=cbind(xobs,yobs)



# 2D random field (function f) ------------------------------
a1=-1.5
a2=0.4

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) - 2
  
}

# exact solution
sol_exact=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact[i] <- z(loc[i,])
}

nnodes = dim(mesh$nodes)[1]
sol_nodes = numeric(nnodes)
for(i in 1:nnodes){
  sol_nodes[i] = z(mesh$nodes[i,])
}


range(sol_exact) 
param = sol_exact
mu<-inv.link(param)
range(mu)
# sampling response:
set.seed(95)
response <- rexp(nloc, rate = 1/mu)


output_CPP <- fdaPDE::smooth.FEM(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = NULL, GCV=GCVFLAG, GCVmethod = GCVmethod,
                                 lambda = lambda, max.steps.FPIRLS=15, family=FAMILY, mu0=NULL, scale.param=NULL)


plot(log10(lambda),output_CPP$GCV)

best_lambda = output_CPP$bestlambda
best_fit = output_CPP$fit.FEM$coeff[,best_lambda]

image(FEM(sol_nodes, FEMbasis))
image(FEM(best_fit, FEMbasis))
