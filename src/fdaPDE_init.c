#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP Density_Estimation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP eval_FEM_fd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP eval_FEM_time(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP eval_FEM_time_nodes(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_mass_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_PDE_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_PDE_space_varying_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_stiff_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_integration_points(SEXP, SEXP, SEXP, SEXP);
extern SEXP points_projection(SEXP, SEXP);
extern SEXP R_triangulate_native(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_Laplace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_Laplace_time(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE_space_varying(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE_space_varying_time(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE_time(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Smooth_FPCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tree_mesh_construction(SEXP, SEXP, SEXP, SEXP);
extern SEXP gam_Laplace(SEXP, SEXP, SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gam_PDE(SEXP, SEXP, SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gam_PDE_space_varying( SEXP, SEXP, SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"Density_Estimation",                (DL_FUNC) &Density_Estimation,                19},
    {"eval_FEM_fd",                       (DL_FUNC) &eval_FEM_fd,                       10},
    {"eval_FEM_time",                     (DL_FUNC) &eval_FEM_time,                     13},
    {"eval_FEM_time_nodes",               (DL_FUNC) &eval_FEM_time_nodes,                5},
    {"get_FEM_mass_matrix",               (DL_FUNC) &get_FEM_mass_matrix,                4},
    {"get_FEM_PDE_matrix",                (DL_FUNC) &get_FEM_PDE_matrix,                17},
    {"get_FEM_PDE_space_varying_matrix",  (DL_FUNC) &get_FEM_PDE_space_varying_matrix,  18},
    {"get_FEM_stiff_matrix",              (DL_FUNC) &get_FEM_stiff_matrix,               4},
    {"get_integration_points",            (DL_FUNC) &get_integration_points,             4},
    {"points_projection",                 (DL_FUNC) &points_projection,                  2},
    {"R_triangulate_native",              (DL_FUNC) &R_triangulate_native,               8},
    {"regression_Laplace",                (DL_FUNC) &regression_Laplace,                18},
    {"regression_Laplace_time",           (DL_FUNC) &regression_Laplace_time,           24},
    {"regression_PDE",                    (DL_FUNC) &regression_PDE,                    21},
    {"regression_PDE_space_varying",      (DL_FUNC) &regression_PDE_space_varying,      22},
    {"regression_PDE_space_varying_time", (DL_FUNC) &regression_PDE_space_varying_time, 28},
    {"regression_PDE_time",               (DL_FUNC) &regression_PDE_time,               27},
    {"Smooth_FPCA",                       (DL_FUNC) &Smooth_FPCA,                       15},
    {"tree_mesh_construction",            (DL_FUNC) &tree_mesh_construction,             4},
    {"gam_Laplace",                       (DL_FUNC) &gam_Laplace,                       24},
    {"gam_PDE",                           (DL_FUNC) &gam_PDE,                           27},
    {"gam_PDE_space_varying",             (DL_FUNC) &gam_PDE_space_varying,             28},
    {NULL, NULL, 0}
};

void R_init_fdaPDE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}