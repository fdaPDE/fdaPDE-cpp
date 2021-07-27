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
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, t_type>::value,
	std::pair<MatrixXr, output_Data<2>> >::type optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, t_type>::value,
	std::pair<MatrixXr, output_Data<2>> >::type optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);

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

	std::pair<MatrixXr, output_Data<2>> solution_bricks;	// Prepare solution to be filled
	
	// Build the Carrier according to problem type
	if(regression.isSV())
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			//Rprintf("Temporal, Areal-forced\n");
			Carrier<InputHandler,Temporal,Forced,Areal>
				carrier = CarrierBuilder<InputHandler>::build_temporal_forced_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal,Forced,Areal>>(carrier);
		}
		else
		{
			//Rprintf("Temporal, Pointwise-forced\n");
			Carrier<InputHandler,Temporal,Forced>
				carrier = CarrierBuilder<InputHandler>::build_temporal_forced_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal,Forced>>(carrier);
		}
	}
	else
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			//Rprintf("Temporal, Areal\n");
			Carrier<InputHandler,Temporal,Areal>
				carrier = CarrierBuilder<InputHandler>::build_temporal_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal,Areal>>(carrier);
		}
		else
		{
			//Rprintf("Temporal, Pointwise\n");
			Carrier<InputHandler,Temporal>
				carrier = CarrierBuilder<InputHandler>::build_temporal_plain_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Temporal>>(carrier);
		}
	}
	
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
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, t_type>::value,
std::pair<MatrixXr, output_Data<2>> >::type optimizer_method_selection(CarrierType & carrier)
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
	else*/ if(optr->get_loss_function() == "GCV" && (optr->get_DOF_evaluation() == "stochastic" || optr->get_DOF_evaluation() == "not_required"))
	{
		//Rprintf("GCV stochastic\n");
		GCV_Stochastic<CarrierType, 2> optim(carrier, true);
		return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 2>, CarrierType>(optim, carrier);
	}
	else // if(optr->get_loss_function() == "unused" && optr->get_DOF_evaluation() == "not_required")
	{
		//Rprintf("Pure evaluation\n");
		GCV_Stochastic<CarrierType, 2> optim(carrier, false);

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// Get the solution
		output_Data<2> output;
		UInt lambdas_count = carrier.get_opt_data()->get_size_S()*carrier.get_opt_data()->get_size_T(); // Valutare se spostare nel Carrier dentro estensione Temporal
		output.z_hat.resize(carrier.get_psip()->rows(), lambdas_count);
		output.lambda_vec.reserve(lambdas_count);
		MatrixXr solution;
 		MatrixXv betas;
		betas.resize(lambdas_count, 1);

		for(UInt j=0; j<carrier.get_opt_data()->get_size_T(); j++)
		{
			for(UInt i=0; i<carrier.get_opt_data()->get_size_S(); i++)
			{
				Real lambdaS = carrier.get_opt_data()->get_lambda_S()[i];
				Real lambdaT = carrier.get_opt_data()->get_lambda_T()[j];
				lambda_type<2> lambda = (lambda_type<2>() << lambdaS, lambdaT).finished();
				output.lambda_vec.push_back(lambda);
				UInt couple_index = j*carrier.get_opt_data()->get_size_S()+i;
				if(i==0 && j==0)
				{
					MatrixXr sol = carrier.apply(lambda); //apply_iterative
					solution.resize(sol.rows(),lambdas_count);
					solution.col(couple_index) = sol;
				}
				else
				{
					solution.col(couple_index) = carrier.apply(lambda); //apply_iterative
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
template<typename EvaluationType, typename CarrierType>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, t_type>::value,
std::pair<MatrixXr, output_Data<2>> >::type optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
{
	// Build wrapper and newton method
	Function_Wrapper<lambda_type<2>, Real, lambda_type<2>, MatrixXr, EvaluationType> Fun(optim);
	typedef Function_Wrapper<lambda_type<2>, Real, lambda_type<2>, MatrixXr, EvaluationType> FunWr;

	const OptimizationData * optr = carrier.get_opt_data();
	//if(optr->get_criterion() == "grid")
	//{
		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// this will be used when grid will be correctly implemented, also for return elements
		
		UInt lambdas_count = carrier.get_opt_data()->get_size_S()*carrier.get_opt_data()->get_size_T();
		std::vector<lambda_type<2>> lambda_vec;
		lambda_vec.reserve(lambdas_count);
		for(UInt j=0; j<carrier.get_opt_data()->get_size_T(); j++)
		{
			Real lambdaT = carrier.get_opt_data()->get_lambda_T()[j];
			for(UInt i=0; i<carrier.get_opt_data()->get_size_S(); i++)
			{
				Real lambdaS = carrier.get_opt_data()->get_lambda_S()[i];
				lambda_type<2> lambda = (lambda_type<2>() << lambdaS, lambdaT).finished();
				lambda_vec.push_back(lambda);
			}
		}
		
		Eval_GCV<lambda_type<2>, MatrixXr, EvaluationType> eval(Fun, lambda_vec);
		output_Data<2> output = eval.Get_optimization_vectorial();

		// Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		MatrixXr solution = carrier.apply(output.lambda_sol); //apply_iterative

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

                //postponed after apply in order to have betas computed
                output.betas = carrier.get_model()->getBeta();

                return {solution, output};

	//}
	//else // 'not_required' optimization can't enter here!! [checked in R code]
	//{
		/*std::unique_ptr<Opt_methods<lambda_type<2>,MatrixXr,EvaluationType>> optim_p =
			Opt_method_factory<lambda_type<2>,MatrixXr,EvaluationType>::create_Opt_method(optr->get_criterion(), Fun);

                // Compute optimal lambda
		Checker ch;
		std::vector<lambda_type<2>> lambda_v_;
		std::vector<Real> GCV_v_;
		Real lambdaS = optr->get_initial_lambda_S();
		Real lambdaT = optr->get_initial_lambda_T();

		if(lambdaS<=0) lambdaS = -1.0;
		if(lambdaT<=0) lambdaT = -1.0;
		
		lambda_type<2> lambda = (lambda_type<2>() << lambdaS, lambdaT).finished();

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		std::pair<lambda_type<2>, UInt> lambda_couple = optim_p->compute(lambda, optr->get_stopping_criterion_tol(), 40, ch, GCV_v_, lambda_v_);

		//Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		// to compute f and g hat
		MatrixXr solution = carrier.apply(lambda_couple.first); //apply_iterative (potrebbe esistere nel GCV un apply_to_b iterative?)

		// postponed after apply in order to have betas computed
		// now the last values in GCV_exact are the correct ones, related to the last iteration
		output_Data<2> output = Fun.get_output(lambda_couple, T, GCV_v_, lambda_v_, ch.which());
		// the copy is necessary for the bulders outside

		return {solution, output};
	}*/
}

#endif
