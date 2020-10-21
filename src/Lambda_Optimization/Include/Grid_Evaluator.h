#ifndef __BATCH_EVALUATOR_H__
#define __BATCH_EVALUATOR_H__

// HEADERS
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"
#include "Function_Variadic.h"
#include "Solution_Builders.h"

// CLASSES
//! Father class for a scalar function evaluation of a given vector of lambda values computing the minimum fuction value
/*!
 \tparam Tuple image type of the gradient of the function
 \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 \tparam Extensions input class if the computations need members already stored in a class
*/
template <typename Tuple, typename Hessian, typename... Extensions>
class Vec_evaluation
{
        protected:
                std::vector<Tuple> lambda_vec;  //!< Vector of lambda to be evaluated

                //! Constructor
                /*!
                \param F_ the function wrapper F performing the evaluation
                \param lambda_vec_ the lambda_vec, vector of lambdas to be evaluated
                */
                Vec_evaluation(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_, const std::vector<Tuple> & lambda_vec_):
                lambda_vec(lambda_vec_), F(F_)
                {
                        // Debugging purpose
                        // Rprintf("Vector evaluator built\n");
                };

                // Function to compute particular parameters related to the mimizing solution.
                /*!
                 It does nothing if not implemented. It is not pure virtual in order to be general
                 and leave the possibility of instantiating the object without implementing that function
                */
                virtual void compute_specific_parameters(void) {};


                //! Only for minimizing solutions. It does nothing if not implemented
                virtual void compute_specific_parameters_best(void) {};

                //! Function evaluator
                Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F;

        public:
                //! Main method function
                /*!
                 \return std::pair<std::vector<Real>, UInt> the vector of evaluations of GCV and the index of the corresponding minimum
                */
                std::pair<std::vector<Real>, UInt> compute_vector(void)
                {
                        UInt dim = lambda_vec.size();
                        UInt index_min = 0; //Assume the first one is the minimum
                        std::vector<Real> evaluations(dim);

                        for (UInt i=0; i<dim; i++)
                        {
                                this->F.set_index(i);
                                evaluations[i] = this->F.evaluate_f(this->lambda_vec[i]); //only scalar functions;

                                this->compute_specific_parameters();
                                if (i==0)
                                     this->compute_specific_parameters_best();

                                if (evaluations[i]<evaluations[index_min])
                                {
                                        this->compute_specific_parameters_best();
                                        index_min=i;
                                }
                        }

                        return {evaluations,index_min};
                }
};

/*!
 Class inheriting form class Vec_evaluation, the function used is GCV evaluation
 \tparam Tuple image type of the gradient of the function
 \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 \tparam Extensions input class if the computations need members already stored in a class
*/
template <typename ...Extensions>
class Eval_GCV: public Vec_evaluation<Real, Real, Extensions...>
{
        protected:
                //! Computes specific parameters needed for GCV
                void compute_specific_parameters(void) override
                {
                        this->F.set_output_partial();
                }

                //! Computes specific parameters needed for GCV best values
                void compute_specific_parameters_best(void) override
                {
                        // Debugging purpose
                        // Rprintf("Specific parameters for GCV computed\n");

                        this->F.set_output_partial_best();
                 }

        public:
                //! Constructor
                /*!
                \param F_ the function wrapper F performing the evaluation
                \param lambda_vec_ the lambda_vec, vector of lambdas to be evaluated
                */
                Eval_GCV(Function_Wrapper<Real, Real, Real, Real, Extensions...> & F_, const std::vector<Real> & lambda_vec_):
                        Vec_evaluation<Real, Real, Extensions...>(F_,lambda_vec_) {};

                //! Function to build the output data
                /*!
                 \return output_Data which contains the almost complete output to be returned to R
                */
                output_Data  Get_optimization_vectorial(void)
                {
                        std::pair<std::vector<Real>, UInt> p = this->compute_vector();
                        output_Data output=this->F.get_output_full();
                        output.GCV_evals  = p.first;
                        output.lambda_sol = this->lambda_vec.at(p.second);      // Safer use of at instead of []
                        output.lambda_pos = 1+p.second;                         // In R numbering
                        output.lambda_vec = this->lambda_vec;
                        output.GCV_opt    = p.first.at(p.second);

                        return output;
                }
};

#endif
