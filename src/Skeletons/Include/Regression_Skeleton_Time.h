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
#include "../../Global_Utilities/Include/Lambda.h"

template<typename CarrierType>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Temporal, CarrierType>::value>, t_type>::value,
	std::pair<MatrixXr, output_Data<2>> >::type optimizer_method_selection(CarrierType & carrier);

template<typename EvaluationType, typename CarrierType, UInt size>
typename std::enable_if<size==2, std::pair<MatrixXr, output_Data<2>>>::type
	optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);
	
template<class GCV_type, typename CarrierType>
std::pair<MatrixXr, output_Data<2>> parabolic_routine(CarrierType & carrier);

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
	
	if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "exact")
	{
		//Rprintf("GCV exact\n");
		if(carrier.get_flagParabolic())	
			return parabolic_routine<GCV_Exact<CarrierType, 1>>(carrier);
		else
		{
			GCV_Exact<CarrierType, 2> optim(carrier);
			return optimizer_strategy_selection<GCV_Exact<CarrierType, 2>, CarrierType, 2>(optim, carrier);
		}

	}
	else if(optr->get_loss_function() == "GCV" && (optr->get_DOF_evaluation() == "stochastic" || optr->get_DOF_evaluation() == "not_required"))
	{		
		//Rprintf("GCV stochastic\n");
		if(carrier.get_flagParabolic())
			return parabolic_routine<GCV_Stochastic<CarrierType, 1>>(carrier);
		else
		{
			GCV_Stochastic<CarrierType, 2> optim(carrier, true);
			return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 2>, CarrierType, 2>(optim, carrier);
		}
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
		output.size_S = carrier.get_opt_data()->get_size_S();
		output.size_T = carrier.get_opt_data()->get_size_T();
		UInt lambdas_count = carrier.get_opt_data()->get_size_S()*carrier.get_opt_data()->get_size_T();
		output.z_hat.resize(carrier.get_psip()->rows(), lambdas_count);
		output.lambda_vec.reserve(lambdas_count);
		output.lambda_vec.clear();
		output.GCV_evals.resize(lambdas_count, -1);
		output.dof.resize(lambdas_count, -1);
		MatrixXr solution;
 		MatrixXv betas;
		betas.resize(lambdas_count, 1);
		for(UInt j=0; j<carrier.get_opt_data()->get_size_T(); j++)
		{
			for(UInt i=0; i<carrier.get_opt_data()->get_size_S(); i++)
			{
				Real lambdaS = carrier.get_opt_data()->get_lambda_S()[i];
				Real lambdaT = carrier.get_opt_data()->get_lambda_T()[j];
				lambda::type<2> lambda = lambda::make_pair(lambdaS, lambdaT);
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

template<class GCV_type, typename CarrierType>
std::pair<MatrixXr, output_Data<2>> parabolic_routine(CarrierType & carrier)
{
	timer Time_partial;
	Time_partial.start();

	std::pair<MatrixXr, output_Data<1>> opt_res;
	std::pair<MatrixXr, output_Data<2>> res;
	
	const OptimizationData * optr = carrier.get_opt_data();

	std::vector<Real> lambdaS = optr->get_lambda_S();
	std::vector<Real> lambdaT = optr->get_lambda_T();
	if(optr->get_criterion() != "grid"){	// set initial data for Newton
		lambdaS[0] = optr->get_initial_lambda_S();
		lambdaT[0] = optr->get_initial_lambda_T();
	}
	res.second.lambda_vec.reserve(res.second.size_S*res.second.size_T);
	res.second.GCV_evals.reserve(res.second.size_S*res.second.size_T);
	res.second.lambda_vec.clear();
	res.second.GCV_evals.clear();
				
	for(UInt j=0; j<optr->get_size_T(); j++)
	{			
		GCV_type optim(carrier, lambdaT[j]);
		std::pair<MatrixXr, output_Data<1>> par_res = 
			optimizer_strategy_selection<GCV_type, CarrierType, 1>(optim, carrier);

		UInt dim = par_res.second.lambda_vec.size();

		for(UInt i=0; i<dim; i++)
			res.second.lambda_vec.push_back(
				lambda::make_pair(par_res.second.lambda_vec[i], lambdaT[j]));
		
		if(par_res.second.GCV_opt < opt_res.second.GCV_opt || j == 0)
		{
			opt_res = par_res;
			res.second.lambda_sol =
				lambda::make_pair(par_res.second.lambda_sol, lambdaT[j]);
			res.second.lambda_pos = j*optr->get_size_S()+par_res.second.lambda_pos;
			res.second.termination = par_res.second.termination;
		}
		
		res.second.rmse.insert(res.second.rmse.end(), par_res.second.rmse.begin(), par_res.second.rmse.end());
		res.second.dof.insert(res.second.dof.end(), par_res.second.dof.begin(), par_res.second.dof.end());
		res.second.GCV_evals.insert(res.second.GCV_evals.end(), par_res.second.GCV_evals.begin(), par_res.second.GCV_evals.end());
		
		res.second.n_it += par_res.second.n_it;
	}

	if(optr->get_criterion() == "grid")
		res.second.size_S = optr->get_size_S();
	else
		res.second.size_S = res.second.lambda_vec.size();  // size_S is the number of lamdaS visited by Newton method
	res.second.size_T = optr->get_size_T();
	res.first = opt_res.first;
	res.second.content = opt_res.second.content;
	timespec T = Time_partial.stop();
	res.second.time_partial = T.tv_sec + 1e-9*T.tv_nsec;
	res.second.z_hat = opt_res.second.z_hat;
	res.second.sigma_hat_sq = opt_res.second.sigma_hat_sq;
	res.second.betas = opt_res.second.betas;
	res.second.GCV_opt = opt_res.second.GCV_opt;
	
	return res;
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
typename std::enable_if<size==2, std::pair<MatrixXr, output_Data<2>>>::type
optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
{
	// Build wrapper and newton method
	Function_Wrapper<lambda::type<2>, Real, lambda::type<2>, MatrixXr, EvaluationType> Fun(optim);
	typedef Function_Wrapper<lambda::type<2>, Real, lambda::type<2>, MatrixXr, EvaluationType> FunWr;
	
	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_criterion() == "grid")
	{
		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		// this will be used when grid will be correctly implemented, also for return elements
		
		UInt lambdas_count = carrier.get_opt_data()->get_size_S()*carrier.get_opt_data()->get_size_T();
		std::vector<lambda::type<2>> lambda_vec;
		lambda_vec.reserve(lambdas_count);
		for(UInt j=0; j<carrier.get_opt_data()->get_size_T(); j++)
		{
			Real lambdaT = carrier.get_opt_data()->get_lambda_T()[j];
			for(UInt i=0; i<carrier.get_opt_data()->get_size_S(); i++)
			{
				Real lambdaS = carrier.get_opt_data()->get_lambda_S()[i];
				lambda::type<2> lambda = lambda::make_pair(lambdaS, lambdaT);
				lambda_vec.push_back(lambda);
			}
		}
		
		Eval_GCV<lambda::type<2>, MatrixXr, EvaluationType> eval(Fun, lambda_vec);
		output_Data<2> output = eval.Get_optimization_vectorial();

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
		std::unique_ptr<Opt_methods<lambda::type<2>,MatrixXr,EvaluationType>> optim_p =
			Opt_method_factory<lambda::type<2>,MatrixXr,EvaluationType>::create_Opt_method(optr->get_criterion(), Fun);			

		// Choose initial lambdaS, lambdaT with grid
		Real lambdaS = optr->get_initial_lambda_S();
		Real lambdaT = optr->get_initial_lambda_T();

		lambda::type<2> lambda_init = lambda::make_pair(lambdaS, lambdaT);
		std::vector<lambda::type<2>> lambda_vec;
		if(lambda_init(0)>0 && lambda_init(1)>0){
			lambda_vec.reserve(5);
			lambda_vec.push_back(lambda_init);
		}
		else
			lambda_vec.reserve(4);
		lambda_vec.push_back(lambda::make_pair(exp(-3), exp(-9)));
		lambda_vec.push_back(lambda::make_pair(exp(-5), exp(-7)));
		lambda_vec.push_back(lambda::make_pair(exp(-7), exp(-5)));
		lambda_vec.push_back(lambda::make_pair(exp(-9), exp(-3)));

		
		UInt dim = lambda_vec.size();
		lambda::type<2> lambda_min;
		Real GCV_min = -1.0;

		for (UInt i=0; i<dim; i++)
		{
		        Rprintf("Pre-Newton grid: evaluating %d/%d\n", i+1, dim);
		        Real evaluation = Fun.evaluate_f(lambda_vec[i]); //only scalar functions;

		        if (evaluation<GCV_min || i==0)
		        {
		      		GCV_min = evaluation;
		      		lambda_min = lambda_vec[i];
		        }
		}

		// Compute optimal lambda
		Checker ch;
		std::vector<lambda::type<2>> lambda_v_;
		std::vector<Real> GCV_v_;

		timer Time_partial; // Of the sole optimization
		Time_partial.start();
		// Rprintf("WARNING: start taking time\n");

		std::pair<lambda::type<2>, UInt> lambda_couple = optim_p->compute(lambda_min, optr->get_stopping_criterion_tol(), 40, ch, GCV_v_, lambda_v_);

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
	}
}

#endif
