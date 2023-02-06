## load all required modules

## regularizing PDES
loadModule("Laplacian_2D_Order1", TRUE)
loadModule("ConstantCoefficients_2D_Order1", TRUE)
loadModule("SpaceVarying_2D_Order1", TRUE)

## SRPDE
loadModule("SRPDE_Laplacian_2D_GeoStatNodes", TRUE)
loadModule("SRPDE_Laplacian_2D_GeoStatLocations", TRUE)
loadModule("SRPDE_Laplacian_2D_Areal", TRUE)
loadModule("SRPDE_ConstantCoefficients_2D_GeoStatNodes", TRUE)
loadModule("SRPDE_ConstantCoefficients_2D_GeoStatLocations", TRUE)
loadModule("SRPDE_ConstantCoefficients_2D_Areal", TRUE)
loadModule("SRPDE_SpaceVarying_2D_GeoStatNodes", TRUE)
loadModule("SRPDE_SpaceVarying_2D_GeoStatLocations", TRUE)
loadModule("SRPDE_SpaceVarying_2D_Areal", TRUE)
