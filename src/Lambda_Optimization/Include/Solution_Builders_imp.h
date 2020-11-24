#ifndef __SOLUTION_BUILDERS_IMP_H__
#define __SOLUTION_BUILDERS_IMP_H__

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP Solution_Builders::build_solution_plain_regression(const MatrixXr & solution, const output_Data & output, const MeshHandler<ORDER, mydim, ndim> & mesh , const InputHandler & regressionData )
{
        // ---- Preparation ----
        // Prepare regresion coefficients space
        MatrixXv beta;
        if(regressionData.getCovariates()->rows()==0)
        {
                beta.resize(1,1);
                beta(0,0).resize(1);
                beta(0,0)(0) = 10e20;
        }
        else
        {
                beta = output.betas;
        }

        // Define string for optimzation method
        UInt code_string;
        if(output.content == "full_optimization")
        {
                code_string = 0;
        }
        else if(output.content == "full_dof_grid")
        {
                code_string = 1;
        }
        else
        {
                code_string = 2;
        }

        const MatrixXr & barycenters = regressionData.getBarycenters();
        const VectorXi & elementIds = regressionData.getElementIds();

        // ---- Copy results in R memory ----
        SEXP result = NILSXP;  // Define emty term --> never pass to R empty or is "R session aborted"
        result = PROTECT(Rf_allocVector(VECSXP, 22)); // 22 elements to be allocated

        // Add solution matrix in position 0
        SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution.rows(), solution.cols()));
        Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < solution.cols(); j++)
	{
		for(UInt i = 0; i < solution.rows(); i++)
			rans[i + solution.rows()*j] = solution(i,j);
	}

        // Add prediction in locations
        SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, output.z_hat.rows(), output.z_hat.cols()));
        rans = REAL(VECTOR_ELT(result, 1));
        for(UInt j = 0; j < output.z_hat.cols(); j++)
        {
                for(UInt i = 0; i < output.z_hat.rows(); i++)
                        rans[i + output.z_hat.rows()*j] = output.z_hat(i,j);
        }

        // Add rmse
        UInt size_rmse = output.rmse.size();
        SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, size_rmse));
        rans = REAL(VECTOR_ELT(result, 2));
        for(UInt j = 0; j < size_rmse; j++)
        {
               rans[j] = output.rmse[j];
        }

        // Add predicted variance of error
        SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, 1));
        rans= REAL(VECTOR_ELT(result, 3));
        rans[0] = output.sigma_hat_sq;

        // Add best lambda value
        SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 4));
        rans[0] = output.lambda_sol;

        // Add best lambda position
        SET_VECTOR_ELT(result, 5, Rf_allocVector(INTSXP, 1));
        UInt *rans1 = INTEGER(VECTOR_ELT(result, 5));
        rans1[0] = output.lambda_pos;

        // Add GCV value
        SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 6));
        rans[0] = output.GCV_opt;

        // Add number of iterations
        SET_VECTOR_ELT(result, 7, Rf_allocVector(INTSXP, 1));
        UInt * rans2 = INTEGER(VECTOR_ELT(result, 7));
        rans2[0] = output.n_it;

        // Add termination criterion
        SET_VECTOR_ELT(result, 8, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 8));
        rans[0] = output.termination;

        // Add code to identify the optimization method
        SET_VECTOR_ELT(result, 9, Rf_allocVector(INTSXP, 1));
        UInt *rans3 = INTEGER(VECTOR_ELT(result, 9));
        rans3[0] = code_string;

        // Add dofs
        UInt size_dof = output.dof.size();
        SET_VECTOR_ELT(result, 10, Rf_allocVector(REALSXP, size_dof));
        rans = REAL(VECTOR_ELT(result, 10));
        for(UInt j = 0; j < size_dof; j++)
        {
               rans[j] = output.dof[j];
        }

        // Add the vector of lambdas
        UInt size_lambda = output.lambda_vec.size();
        SET_VECTOR_ELT(result, 11, Rf_allocVector(REALSXP, size_lambda));
        rans = REAL(VECTOR_ELT(result, 11));
        for(UInt j = 0; j < size_lambda; j++)
        {
               rans[j] = output.lambda_vec[j];
        }

        // Add vector of GCV
        UInt size_vec = output.GCV_evals.size();
        SET_VECTOR_ELT(result, 12, Rf_allocVector(REALSXP, size_vec));
        rans = REAL(VECTOR_ELT(result, 12));
        for(UInt j = 0; j < size_vec; j++)
        {
               rans[j] = output.GCV_evals[j];
        }

        // Add time employed
        SET_VECTOR_ELT(result, 13, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 13));
        rans[0] = output.time_partial;

        // Copy betas
        SET_VECTOR_ELT(result, 14, Rf_allocMatrix(REALSXP, beta(0).size(), beta.size()));
        Real *rans4 = REAL(VECTOR_ELT(result, 14));
        for(UInt j = 0; j < beta.size(); j++)
        {
                for(UInt i = 0; i < beta(0).size(); i++)
                        rans4[i + beta(0).size()*j] = beta(j)(i);
        }

        if(regressionData.getSearch()==2){
            // Send tree information to R
            SET_VECTOR_ELT(result, 15, Rf_allocVector(INTSXP, 1)); //tree_header information
            int *rans5 = INTEGER(VECTOR_ELT(result, 15));
            rans5[0] = mesh.getTree().gettreeheader().gettreelev();

            SET_VECTOR_ELT(result, 16, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
            Real *rans6 = REAL(VECTOR_ELT(result, 16));
            for(UInt i = 0; i < ndim*2; i++)
            rans6[i] = mesh.getTree().gettreeheader().domainorig(i);

            SET_VECTOR_ELT(result, 17, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
            Real *rans7 = REAL(VECTOR_ELT(result, 17));
            for(UInt i = 0; i < ndim*2; i++)
            rans7[i] = mesh.getTree().gettreeheader().domainscal(i);


            UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
            SET_VECTOR_ELT(result, 18, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
            int *rans8 = INTEGER(VECTOR_ELT(result, 18));
            for(UInt i = 0; i < num_tree_nodes; i++)
                    rans8[i] = mesh.getTree().gettreenode(i).getid();

            for(UInt i = 0; i < num_tree_nodes; i++)
                    rans8[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

            for(UInt i = 0; i < num_tree_nodes; i++)
                    rans8[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

            SET_VECTOR_ELT(result, 19, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
            Real *rans9 = REAL(VECTOR_ELT(result, 19));
            for(UInt j = 0; j < ndim*2; j++)
            {
                    for(UInt i = 0; i < num_tree_nodes; i++)
                            rans9[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
            }
        }

        // Send barycenter information to R
        SET_VECTOR_ELT(result, 20, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
        int *rans10 = INTEGER(VECTOR_ELT(result, 20));
        for(UInt i = 0; i < elementIds.rows(); i++)
        rans10[i] = elementIds(i);

        SET_VECTOR_ELT(result, 21, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
        Real *rans11 = REAL(VECTOR_ELT(result, 21));
        for(UInt j = 0; j < barycenters.cols(); j++)
        {
                for(UInt i = 0; i < barycenters.rows(); i++)
                        rans11[i + barycenters.rows()*j] = barycenters(i,j);
        }

        UNPROTECT(1);

        return(result);
}

#endif
