#ifndef __REGRESSION_SKELETON_H__
#define __REGRESSION_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../Lambda_Optimization/Include/Carrier.h"
#include "../../Lambda_Optimization/Include/Grid_Evaluator.h"
#include "../../Lambda_Optimization/Include/Lambda_Optimizer.h"
#include "../../Lambda_Optimization/Include/Newton.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Methods_Factory.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"

template<typename CarrierType>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, f_type>::value,
	std::pair<MatrixXr, output_Data<1>> >::type optimizer_method_selection(CarrierType & carrier);

template<typename EvaluationType, typename CarrierType, UInt size>
typename std::enable_if<size==1, std::pair<MatrixXr, output_Data<1>>>::type 
	optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler & regressionData, OptimizationData & optimizationData, SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, regressionData.getSearch());	// Create the mesh
	MixedFERegression<InputHandler> regression(regressionData, optimizationData, mesh.num_nodes()); // Define the mixed object

	regression.preapply(mesh); // preliminary apply (preapply) to store all problem matrices

    std::pair<MatrixXr, output_Data<1>> solution_bricks;	// Prepare solution to be filled

	// Build the Carrier according to problem type
	if(regression.isSV())
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			//Rprintf("Areal-forced\n");
			Carrier<InputHandler,Forced,Areal>
				carrier = CarrierBuilder<InputHandler>::build_forced_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler, Forced,Areal>>(carrier);
		}
		else
		{
			//Rprintf("Pointwise-forced\n");
			Carrier<InputHandler,Forced>
				carrier = CarrierBuilder<InputHandler>::build_forced_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Forced>>(carrier);
		}
	}
	else
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			//Rprintf("Areal\n");
			Carrier<InputHandler,Areal>
				carrier = CarrierBuilder<InputHandler>::build_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Areal>>(carrier);
		}
		else
		{
			//Rprintf("Pointwise\n");
			Carrier<InputHandler>
				carrier = CarrierBuilder<InputHandler>::build_plain_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler>>(carrier);
		}
	}

 	return Solution_Builders::build_solution_plain_regression<InputHandler, ORDER, mydim, ndim>(solution_bricks.first, solution_bricks.second, mesh, regressionData, regression);
}

//! Function to select the right optimization method
/*
 \tparam CarrierType the type of Carrier to be employed
 \param carrier the Carrier used for the methods
 \return the solution to pass to the Solution_Builders
*/
template<typename CarrierType>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, f_type>::value,
std::pair<MatrixXr, output_Data<1>> >::type optimizer_method_selection(CarrierType & carrier)
{
	// Build the optimizer
	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "exact")
	{
		//Rprintf("GCV exact\n");
		GCV_Exact<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Exact<CarrierType, 1>, CarrierType, 1>(optim, carrier);
	}
	else if(optr->get_loss_function() == "GCV" && (optr->get_DOF_evaluation() == "stochastic" || optr->get_DOF_evaluation() == "not_required"))
	{
		//Rprintf("GCV stochastic\n");
		GCV_Stochastic<CarrierType, 1> optim(carrier, true);
		return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 1>, CarrierType, 1>(optim, carrier);
	}
	else // if(optr->get_loss_function() == "unused" && optr->get_DOF_evaluation() == "not_required")
	{
		//Rprintf("Pure evaluation\n");
		GCV_Stochastic<CarrierType, 1> optim(carrier, false);

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// Get the solution
		output_Data<1> output;
		output.z_hat.resize(carrier.get_psip()->rows(),carrier.get_opt_data()->get_size_S());
		output.lambda_vec = carrier.get_opt_data()->get_lambda_S();
		MatrixXr solution;
 		MatrixXv betas;
		betas.resize(carrier.get_opt_data()->get_size_S(),1);

		for(UInt j=0; j<carrier.get_opt_data()->get_size_S(); j++)
		{
			if(j==0)
			{
				MatrixXr sol = carrier.apply(carrier.get_opt_data()->get_lambda_S()[j]);
				solution.resize(sol.rows(),carrier.get_opt_data()->get_size_S());
				solution.col(j) = sol;
			}
			else
			{
				solution.col(j) = carrier.apply(carrier.get_opt_data()->get_lambda_S()[j]);
			}
			optim.combine_output_prediction(solution.topRows(solution.rows()/2).col(j),output,j);
			if(carrier.get_model()->getBeta().cols()>0 && carrier.get_model()->getBeta().rows()>0)
				betas.coeffRef(j,0)=carrier.get_model()->getBeta().coeffRef(0,0);
		}

		// Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

                // postponed after apply in order to have betas computed
        output.betas = betas;

        return {solution, output};
	}
}

//! Function to apply the optimization strategy, grid or Newton
/*
\\tparam EvaluationType optimization type to be used
 \tparam CarrierType the type of Carrier to be employed
 \param optim EvaluationType containing the class related to the function to be optimized,together with the method (exact or stochastic)
 \param carrier the Carrier used for the methods
 \return the solution to pass to the Solution_Builders
*/

template<typename EvaluationType, typename CarrierType, UInt size>
typename std::enable_if<size==1, std::pair<MatrixXr, output_Data<1>>>::type 
	optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
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
		output_Data<1> output = eval.Get_optimization_vectorial();

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

        // Choose initial lambdaS with grid
		Real lambdaS_init = optr->get_initial_lambda_S();   // first value of lambdaS sequence

		std::vector<Real> lambdaS_grid = {5.000000e-05, 1.442700e-03, 4.162766e-02, 1.201124e+00, 3.465724e+01, 1.000000e+03};
			// Start from 6 lambda and find the minimum value of GCV to start from it the newton's method


		UInt dim = lambdaS_grid.size();
	        Real lambdaS_min;
	        Real GCV_min = -1.0;

	        for (UInt i=0; i<dim; i++)
	        {
	                Rprintf("Pre-Newton grid: evaluating %d/%d\n", i+1, dim);
	                Real evaluation = Fun.evaluate_f(lambdaS_grid[i]); //only scalar functions;

	                if (evaluation<GCV_min || i==0)
	                {
	              		GCV_min = evaluation;
	              		lambdaS_min = lambdaS_grid[i];
	                }
	        }

		// If lambdaS_init <= 0, use the one from grid
		if (lambdaS_init>lambdaS_min/4 || lambdaS_init<=0)
			lambdaS_init = lambdaS_min/8;

		Checker ch;
		std::vector<Real> lambda_v_;
		std::vector<Real> GCV_v_;
		

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// Compute optimal lambda
		std::pair<Real, UInt> lambda_couple = optim_p->compute(lambdaS_init, optr->get_stopping_criterion_tol(), 40, ch, GCV_v_, lambda_v_);

		//Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		// to compute f and g hat
		MatrixXr solution = carrier.apply(lambda_couple.first);

		// postponed after apply in order to have betas computed
		// now the last values in GCV_exact are the correct ones, related to the last iteration
		output_Data<1> output = Fun.get_output(lambda_couple, T, GCV_v_, lambda_v_, ch.which());
		// the copy is necessary for the bulders outside

		return {solution, output};
	}
}

#endif
