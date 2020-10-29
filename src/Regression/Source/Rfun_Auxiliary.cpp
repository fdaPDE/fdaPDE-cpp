#include "../../FdaPDE.h"
#include "../../Skeletons/Include/Auxiliary_Skeleton.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"

extern "C"
{
        //! A function required for anysotropic and nonstationary regression (only 2D)
        /*!
            \return points where the PDE space-varying params are evaluated in the R code
        */
        SEXP get_integration_points(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
        {
        	//Declare pointer to access data from C++
        	int order = INTEGER(Rorder)[0];

        	//Get mydim and ndim
        	UInt mydim=INTEGER(Rmydim)[0];
        	UInt ndim=INTEGER(Rndim)[0];

            if(order == 1 && ndim==2)
                return(get_integration_points_skeleton<1,2,2>(Rmesh));
            else if(order == 2 && ndim==2)
                return(get_integration_points_skeleton<2,2,2>(Rmesh));
            if(order == 1 && ndim==3 && mydim==2)
                return(get_integration_points_skeleton<1,2,3>(Rmesh));
            else if(order == 2 && ndim==3 && mydim==2)
                return(get_integration_points_skeleton<2,2,3>(Rmesh));
            if(order == 1 && mydim==3)
                return(get_integration_points_skeleton<1,3,3>(Rmesh));
            else if(order == 2 && mydim==3)
                return(get_integration_points_skeleton<2,3,3>(Rmesh));

            return(NILSXP);
        }

        //! A utility, not used for system solution, may be used for debugging
        SEXP get_FEM_mass_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
        {
        	int order = INTEGER(Rorder)[0];

        	//Get mydim and ndim
        	UInt mydim=INTEGER(Rmydim)[0];
        	UInt ndim=INTEGER(Rndim)[0];

        	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);

                if(order==1 && ndim==2)
                        return(get_FEM_Matrix_skeleton<1,2,2>(Rmesh, mass));
                if(order==2 && ndim==2)
                        return(get_FEM_Matrix_skeleton<2,2,2>(Rmesh, mass));
                return(NILSXP);
        }

        //! A utility, not used for system solution, may be used for debugging
        SEXP get_FEM_stiff_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
        {
        	int order = INTEGER(Rorder)[0];

        	//Get mydim and ndim
        	UInt mydim=INTEGER(Rmydim)[0];
        	UInt ndim=INTEGER(Rndim)[0];

        	typedef EOExpr<Stiff> ETMass;   Stiff EStiff;   ETMass stiff(EStiff);

                if(order==1 && ndim==2)
                        return(get_FEM_Matrix_skeleton<1,2,2>(Rmesh, stiff));
                if(order==2 && ndim==2)
                        return(get_FEM_Matrix_skeleton<2,2,2>(Rmesh, stiff));
                return(NILSXP);
        }
}
