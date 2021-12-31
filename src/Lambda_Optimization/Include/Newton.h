#ifndef __NEWTON_H__
#define __NEWTON_H__

// HEADERS
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include "Function_Variadic.h"
#include "../../Global_Utilities/Include/Lambda.h"

// CLASSES
//! Checker, contains data regarding the last process, used in optimization processes
class Checker
{
        private:
                bool reached_max_iter;          //!< Boolean for maximum number ot iterations reached
                bool reached_tolerance;         //!< Boolean for tolerance reached

        public:
                //! Basic Constructor: everything set as false
                Checker(void): reached_max_iter(false), reached_tolerance(false) {}

                //! Sets max number of iterations as true
                inline void set_max_iter(void)  {reached_max_iter  = true;}

                //! Sets the tolerance for the optimization method as true
                inline void set_tolerance(void) {reached_tolerance = true;}

                //! Returns the reason of conclusion of the iterative method
                /*!
                 \return the code type of the problem (1 tolerance, 2 max iterations, -1 error)
                */
                inline int which(void) const
                {
                        if (reached_tolerance == true)
                                return 1;
                        else if (reached_max_iter ==  true)

                                return 2;
                        else
                                return -1; //error status
                }
};

//! Father class to apply generic optimization methods
/*!
 * \tparam Tuple image type of the gradient of the function
 * \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam Extensions input class if the computations need members already stored in a class
 */
template <typename Tuple, typename Hessian, typename... Extensions>
class Opt_methods
{
        protected:
                //! Constructor
                /*!
                 \param F_ the function wrapper F to be optimized
                */
                Opt_methods(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_): F(F_) {}

               //! Function evaluator
                Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F;

        public:

                //! Function to apply the optimization method and obtain as a result the couple (optimal lambda, optimal value of the function)
                virtual std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<Tuple> & lambda_v) = 0;

                //! Virtual Destuctor
                virtual ~Opt_methods(){};
};

template <typename Tuple>
struct Auxiliary
{
        // NOT yet implemented
};

//! Auxiliary class to perform elementary mathematical operations and checks: specialization for 1 dimensional case
template<>
struct Auxiliary<lambda::type<1>>
{
        public:
                Auxiliary(void) {};

                static inline bool isNull(Real n)                       {return (n == 0);}      //!< Check if the input value is zero \param n number to be checked
                static inline void divide(Real a, Real b, Real & x)     {x = b/a;}              //!< Apply a division \param a denominator \param b numerator \param x reference to result
                static inline Real residual(Real a)                     {return std::abs(a);}   //!< Compute the norm of the residual \param a take the absolute value of this number \return the absolute value of a
                static inline bool isPositive(Real x)                   {return x>0;}           //!< Check if the input is positive
};

//! Auxiliary class to perform elementary mathematical operations and checks: specialization for n dimensional case
template<>
struct Auxiliary<lambda::type<2>>
{
        public:
                Auxiliary(void) {};

                //! Check if the input value is zero
                /*!
                 \param n matrix to be checked
                */
                static inline bool isNull(MatrixXr n)
                {
                        UInt sz = n.size();
                        return (n == MatrixXr::Zero(sz,sz));
                }

                //! Solve a linear system in the optimization method
                /*!
                 \param A system matrix
                 \param b right-hand side
                 \param x reference to result vector
                */
                static inline void divide(const MatrixXr & A, const lambda::type<2> & b, lambda::type<2> & x)
                {
                        Cholesky::solve(A, b, x);
                }

                //! Compute the norm of the residual
                /*!
                 \param a vector of which compute the norm
                 \return the computed norm
                */
                static inline Real residual(lambda::type<2> a)
                {
                        return a.norm();
                }

                //! Check if the input is positive
                /*!
                 \param x vector of which to check positivity
                 \return positivity of x
                */
                static inline bool isPositive(lambda::type<2> x)
                {
                        return x(0)>0 && x(1)>0;
                }

};


//! Class to apply Newton exact method, inheriting from Opt_methods
/*!
 * \tparam Tuple image type of the gradient of the function
 * \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam Extensions input class if the computations need members already stored in a class
*/
template <typename Tuple, typename Hessian, typename ...Extensions>
class Newton_ex: public Opt_methods<Tuple, Hessian, Extensions...>
{

};

template <typename ...Extensions>
class Newton_ex<lambda::type<1>, Real, Extensions...>: public Opt_methods<lambda::type<1>, Real, Extensions...>
{
        public:

                // Constructor
                /*!
                 \param F_ the function wrapper F to be optimized
                 \note F cannot be const, it must be modified
                */
                Newton_ex(Function_Wrapper<lambda::type<1>, Real, lambda::type<1>, Real, Extensions...> & F_): Opt_methods<lambda::type<1>, Real, Extensions...>(F_) {};

                //! Apply Newton fd method
                std::pair<Real, UInt> compute(const lambda::type<1> & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<lambda::type<1>> & lambda_v) override;


                //! Virtual Destuctor
                virtual ~Newton_ex(){};


};

template <typename ...Extensions>
class Newton_ex<lambda::type<2>, MatrixXr, Extensions...>: public Opt_methods<lambda::type<2>, MatrixXr, Extensions...>
{
        public:

                // Constructor
                /*!
                 \param F_ the function wrapper F to be optimized
                 \note F cannot be const, it must be modified
                */
                Newton_ex(Function_Wrapper<lambda::type<2>, Real, lambda::type<2>, MatrixXr, Extensions...> & F_): Opt_methods<lambda::type<2>, MatrixXr, Extensions...>(F_) {};

                //! Apply Newton fd method
                std::pair<lambda::type<2>, UInt> compute(const lambda::type<2> & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<lambda::type<2>> & lambda_v) override;


                //! Virtual Destuctor
                virtual ~Newton_ex(){};


};


template <typename Tuple, typename Hessian, typename ...Extensions>
class Newton_fd: public Opt_methods<Tuple, Hessian, Extensions...>
{
        // NOT yet implemented
};

//! Class to apply Newton method exploiting finite differences to compute derivatives, inheriting from Opt_methods
/*!
 \tparam Tuple image type of the gradient of the function
 \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 \tparam Extensions input class if the computations need members already stored in a class
*/
template <typename ...Extensions>
class Newton_fd<lambda::type<1>, Real, Extensions...>: public Opt_methods<lambda::type<1>, Real, Extensions...>
{
        public:

                // Constructor
                /*!
                 \param F_ the function wrapper F to be optimized
                 \note F cannot be const, it must be modified
                */
                Newton_fd(Function_Wrapper<lambda::type<1>, Real, lambda::type<1>, Real, Extensions...> & F_): Opt_methods<lambda::type<1>, Real, Extensions...>(F_) {};

                //! Apply Newton fd method
                std::pair<Real, UInt> compute(const lambda::type<1> & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<lambda::type<1>> & lambda_v) override;


                //! Virtual Destuctor
                virtual ~Newton_fd(){};


};

template <typename ...Extensions>
class Newton_fd<lambda::type<2>, MatrixXr, Extensions...>: public Opt_methods<lambda::type<2>, MatrixXr, Extensions...>
{
        public:

                // Constructor
                /*!
                 \param F_ the function wrapper F to be optimized
                 \note F cannot be const, it must be modified
                */
                Newton_fd(Function_Wrapper<lambda::type<2>, Real, lambda::type<2>, MatrixXr, Extensions...> & F_): Opt_methods<lambda::type<2>, MatrixXr, Extensions...>(F_) {};

                //! Apply Newton fd method
                std::pair<lambda::type<2>, UInt> compute(const lambda::type<2> & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<lambda::type<2>> & lambda_v) override;


                //! Virtual Destuctor
                virtual ~Newton_fd(){};


};


#include "Newton_imp.h"

#endif
