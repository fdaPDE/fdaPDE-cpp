#ifndef __REGRESSION_SKELETON_TIME_H__
#define __REGRESSION_SKELETON_TIME_H__

#include "../../FdaPDE.h"
#include "../../Lambda_Optimization/Include/Carrier.h"
#include "../../Lambda_Optimization/Include/Grid_Evaluator.h"
#include "../../Lambda_Optimization/Include/Lambda_Optimizer.h"
#include "../../Lambda_Optimization/Include/Newton.h"
#include "../../Lambda_Optimization/Include/Optimization_Methods_Factory.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"

template<typename CarrierType>
std::pair<MatrixXr, output_Data<std::pair<Real, Real>>> optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
std::pair<MatrixXr, output_Data<std::pair<Real, Real>>> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton_time(InputHandler & regressionData, OptimizationData & optimizationData, SEXP Rmesh, SEXP Rmesh_time)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, regressionData.getSearch());//! load the mesh
	UInt n_time = Rf_length(Rmesh_time);
	std::vector<Real> mesh_time(n_time);
	for(UInt i=0; i<n_time; ++i)
	{
		mesh_time[i] = REAL(Rmesh_time)[i];
	}
	MixedFERegression<InputHandler> regression(mesh_time, regressionData, optimizationData, mesh.num_nodes());//! load data in a C++ object

	regression.preapply(mesh); // preliminary apply (preapply) to store all problem matrices

    std::pair<MatrixXr, output_Data<std::pair<Real, Real>>> solution_bricks;	// Prepare solution to be filled
/*
    if(regressionData.getFlagParabolic()){
    	//DO forced areal

    }
    else{
    	//DO forced areal SEPARABLE
    }
*/
    Carrier<InputHandler,Parabolic>
				carrier = CarrierBuilder<InputHandler>::build_temporal_carrier(regressionData, regression, optimizationData);
	solution_bricks = optimizer_method_selection<Carrier<InputHandler, Parabolic>>(carrier);
	//Se si riesce, usare plain anche nel caso temporal
 	return Solution_Builders::build_solution_temporal_regression<InputHandler, ORDER, mydim, ndim>(solution_bricks.first, solution_bricks.second, mesh, regressionData, regression);
}



//! Function to select the right optimization method
/*
 \tparam CarrierType the type of Carrier to be employed
 \param carrier the Carrier used for the methods
 \return the solution to pass to the Solution_Builders
*/
template<typename CarrierType>
std::pair<MatrixXr, output_Data<std::pair<Real, Real>>> optimizer_method_selection(CarrierType & carrier)
{
	// Build the optimizer
	const OptimizationData * optr = carrier.get_opt_data();
	/*
	if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "exact")
	{
		//Rprintf("GCV exact\n");
		GCV_Exact<CarrierType, 2> optim(carrier);
		return optimizer_strategy_selection<GCV_Exact<CarrierType, 2>, CarrierType>(optim, carrier);
	}
	else if(optr->get_loss_function() == "GCV" && (optr->get_DOF_evaluation() == "stochastic" || optr->get_DOF_evaluation() == "not_required"))
	{
		//Rprintf("GCV stochastic\n");
		GCV_Stochastic<CarrierType, 2> optim(carrier, true);
		return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 2>, CarrierType>(optim, carrier);
	}
	else // if(optr->get_loss_function() == "unused" && optr->get_DOF_evaluation() == "not_required")
	{
		*/
		//Rprintf("Pure evaluation\n");
		GCV_Stochastic<CarrierType, 2> optim(carrier, false);

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// Get the solution
		output_Data<std::pair<Real, Real>> output;
		UInt lambdas_count = carrier.get_opt_data()->get_size_S()*carrier.get_opt_data()->get_size_T(); // Valutare se spostare nel Carrier dentro estensione Temporal
		output.z_hat.resize(carrier.get_psip()->rows(), lambdas_count);
		output.lambda_vec.reserve(lambdas_count);
		MatrixXr solution;
 		MatrixXv betas;
		betas.resize(lambdas_count, 1);

		for(UInt i=0; i<carrier.get_opt_data()->get_size_S(); i++)
		{
			for(UInt j=0; j<carrier.get_opt_data()->get_size_T(); j++)
			{
				Real lambdaS = carrier.get_opt_data()->get_lambda_S()[i];
				Real lambdaT = carrier.get_opt_data()->get_lambda_T()[j];
				output.lambda_vec.push_back(std::make_pair(lambdaS, lambdaT));
				UInt couple_index = i*carrier.get_opt_data()->get_size_T()+j;
				if(i==0 && j==0)
				{
					MatrixXr sol = carrier.apply(lambdaS, lambdaT);
					solution.resize(sol.rows(),lambdas_count);
					solution.col(couple_index) = sol;
				}
				else
				{
					solution.col(couple_index) = carrier.apply(lambdaS, lambdaT);
				}
				optim.combine_output_prediction(solution.topRows(solution.rows()/2).col(couple_index),output,couple_index);
				if(carrier.get_model()->getBeta().cols()>0 && carrier.get_model()->getBeta().rows()>0)
					betas.coeffRef(couple_index,0) = carrier.get_model()->getBeta().coeffRef(0,0);
			}
		}

		// Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

        // postponed after apply in order to have betas computed
        output.betas = betas;

        return {solution, output};
	//}
}

//! Function to apply the optimization strategy, grid or Newton
/*
\\tparam EvaluationType optimization type to be used
 \tparam CarrierType the type of Carrier to be employed
 \param optim EvaluationType containing the class related to the function to be optimized,together with the method (exact or stochastic)
 \param carrier the Carrier used for the methods
 \return the solution to pass to the Solution_Builders
*/
/*
template<typename EvaluationType, typename CarrierType>
std::pair<MatrixXr, output_Data<Real>> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
{
	// Build wrapper and newton method
	Function_Wrapper<Real, Real, Real, Real, EvaluationType> Fun(optim);
	typedef Function_Wrapper<Real, Real, Real, Real, EvaluationType> FunWr;

	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_criterion() == "grid")
	{
		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// this will be used when grid will be correctly implemented, also for return elements

		Eval_GCV<Real, Real, EvaluationType> eval(Fun, optr->get_lambda_S());
		output_Data<Real> output = eval.Get_optimization_vectorial();

		// Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		MatrixXr solution = carrier.apply(output.lambda_sol);

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

                //postponed after apply in order to have betas computed
                output.betas = carrier.get_model()->getBeta();

                return {solution, output};

	}
	else // 'not_required' optimization can't enter here!! [checked in R code]
	{
		std::unique_ptr<Opt_methods<Real,Real,EvaluationType>> optim_p =
			Opt_method_factory<Real, Real, EvaluationType>::create_Opt_method(optr->get_criterion(), Fun);

                // Compute optimal lambda
		Checker ch;
		std::vector<Real> lambda_v_;
		std::vector<Real> GCV_v_;
		Real lambda = optr->get_initial_lambda_S();

		if(lambda<=0)
		{
			lambda = -1.0;
		}

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		std::pair<Real, UInt> lambda_couple = optim_p->compute(lambda, optr->get_stopping_criterion_tol(), 40, ch, GCV_v_, lambda_v_);

		//Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		// to compute f and g hat
		MatrixXr solution = carrier.apply(lambda_couple.first);

		// postponed after apply in order to have betas computed
		// now the last values in GCV_exact are the correct ones, related to the last iteration
		output_Data<Real> output = Fun.get_output(lambda_couple, T, GCV_v_, lambda_v_, ch.which());
		// the copy is necessary for the bulders outside

		return {solution, output};
	}
}
*/

#endif
