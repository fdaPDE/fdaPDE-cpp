#ifndef __GOF_UPDATER_H__
#define __GOF_UPDATER_H__

// HEADERS
#include <functional>
#include "../../FdaPDE.h"

// CLASSES
//! External updater for Lambda Optimizer
/*!
 This class works as an external updater for a LambdaOptimizer method,
 storing its pointer and keeping track of the updates performed in the past.
 The purpose of this structure is mainly keeping at minimum the complexity of
 each update by performing just what is needed from time to time.
 \tparam LambdaOptim LambdaOptimizer type on which the update has to be performed
 \tparam T type of the optimization parameter lambda
*/
template <typename LambdaOptim, typename T>
class GOF_updater
{
        private:
                // -- PARAMETERS FOR THE UPDATES --
                //! std::vector of lambdas to keep track of the last lambdas used [position means degree of derivative]
                std::vector<T> last_lambda_derivatives;
                //! std::vector storing the updaters to be called [position means degree of derivative]
                std::vector<std::function<void(Real)>> updaters;
                //! pointer collecting the lambda optimizer on which the updates have to be performed
                LambdaOptim * start_ptr = nullptr;

                // -- PRIVATE MEMBERS --
                //! Function that selectively updates just the essential terms for the needed task
                /*!
                 This functions calls all the updaters on the given pointer from start to finish
                 [e.g. start==1, finish==2 means: call first and second updater to prepare
                 Lambda Optimizer for first and second derivative computation]
                 \param start first updater to be called, increase sequentially from this
                 \param finish last updater to be called (included)
                 \param lambda the actual value of lambda to be used for the computation
                */
                inline void call_from_to(UInt start, UInt finish, T lambda)
                {
                        for(UInt i=start; i<=finish; ++i) //loop from start to finish included
                        {
                                updaters[i](lambda);                 // call the update
                                last_lambda_derivatives[i] = lambda; // keep track of the performed update
                        }
                }

                //! Function that initializes the updaters given the pinter to the optimizatio method
                /*!
                 \param lopt_ptr pointer to the Lambda Optimizer from which to capture the updaters
                */
                inline void updaters_setter(LambdaOptim * lopt_ptr)
                {
                        this->updaters.reserve(3); // all methods need up to the second derivative
                        this->updaters.push_back(std::bind(&LambdaOptim::zero_updater, lopt_ptr, std::placeholders::_1));
                        this->updaters.push_back(std::bind(&LambdaOptim::first_updater, lopt_ptr, std::placeholders::_1));
                        this->updaters.push_back(std::bind(&LambdaOptim::second_updater, lopt_ptr, std::placeholders::_1));
                }

        public:
                // -- CONSTRUCTORS --
                //! Default constructor
                GOF_updater(void) = default;

                //! Initializier for the first lambdas to be stored
                /*!
                 \param first_lambdas the first lambdas to be stored in the vector of last_lambda_derivatives
                */
                inline void initialize(const std::vector<T> & first_lambdas)
                {
                        last_lambda_derivatives = first_lambdas;
                }

                // -- PUBLIC UPDATER --
                //! Public function that selectively updates just the essential terms for the needed task
                /*!
                 This functions calls all the updaters on the given pointer up to finish
                 [e.g. finish==2 calls, if not already updated, zero first and second updater]
                 \param finish last updater to be called (included)
                 \param lambda the actual value of lambda to be used for the computation
                 \param lopt_ptr the object from which to perform the update
                */
                inline void call_to(UInt finish, T lambda, LambdaOptim * lopt_ptr)
                {
                        if(start_ptr != lopt_ptr) // new pointer to be stored (or first time we store)
                        {
                                //Debugging purpose
                                //Rprintf("--- Set updaters ---\n");
                                initialize(std::vector<Real>{-1.,-1.,-1.});     // dummy initialize the last lambdas
                                updaters_setter(lopt_ptr);                      // set all the updaters from the given pointer
                                start_ptr = lopt_ptr;                           // keep track of the pointer to avoid this procedure next time
                        }

                        bool found = false;                                     // cycle breaker
                        for(UInt i = 0; i<=finish && found==false; ++i)         // loop until the desired update level (finish)
                                if(lambda != last_lambda_derivatives[i])        // if we have found the first not updated derivative
                                {
                                        call_from_to(i, finish, lambda);        // update from that derivative to the needed level (finish)
                                        found = true;                           // break the cycle since the update is complete
                                }
                }
};

#endif
