#ifndef __FUNCTION_VARIADIC_H__
#define __FUNCTION_VARIADIC_H__

// HEADERS
#include <functional>
#include <iostream>
#include <type_traits>
#include "../../FdaPDE.h"

// CLASSES

//! Function wrapper. Stores a function from R^m to R^n. Variadic template to store either from inheritance or directly.
/*!
This class implements a general function, giving the possibility of computing also the first and second derivatives. It is used to simplify the code to instantiate it, since by Variadic Templates it has the same instantiation whichever class it derives from. The class from which possibly inherits need to have the members compute_f, compute_fp, compute_fs, which are wrapped in this class.
 \tparam Dtype domain type of the function
 \tparam Ctype image type of the function
 \tparam Tuple image type of the gradient of the function
 \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 \tparam Extensions... input class if the computations need members already stored in a class
*/
 template<typename Dtype, typename Ctype, typename Tuple, typename Hessian, typename... Extensions>
 class Function_Wrapper : public Extensions...
 {
        private:
                // Actual functions (function, first derivative/gradient, second derivative/Hessian) stored inside the Function class
                std::function<Ctype(Dtype)>   g;   //!< base function
                std::function<Tuple(Dtype)>   dg;  //!< first derivative
                std::function<Hessian(Dtype)> ddg; //!< second_derivative
                
        public:
                //  -- CONSTRUCTORS --
                //! Constructor taking object Extensions and initializing with it the new class
                template<typename D=typename std::enable_if<sizeof...(Extensions)!=0,void>, typename ...T>
                Function_Wrapper(T&& ...ext):Extensions(std::forward<T>(ext))...{};

                //! Since I have defined a constructor I need to indicate the default constructor
                Function_Wrapper() = default;

                // -- SETTERS --
                //! Function to inizialize the g,dg,ddg, used only if you explicit directly the functions
                /*!
                 \param g_ the std:::function from which to copy
                 \param dg_ the derivative std:::function from which to copy
                 \param ddg_ the second derivative std:::function from which to copy
                */
                inline void set_functions(const std::function<Ctype(Dtype)> & g_,const std::function<Tuple(Dtype)> & dg_,const std::function<Hessian(Dtype)> & ddg_)
                {
                        g   = g_;
                        dg  = dg_;
                        ddg = ddg_;
                }

                // -- OPERATORS --
                //! Function version of a std::function () operator, use SFINAE to differentiate from the two possibilities of instantiation, derived case
                /*!
                 \tparam U used for SFINAE purpose, is the type of lambda parameter, should not be void
                 \param lambda the independent variable for function evaluation, usually a Real or a VectorXr
                 \return a Ctype: evaluation of the function at point lambda
                 \note never choose U void (it is unresonable and beaks the mechanism)
                */
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)==0 || std::is_void<U>::value, Ctype>::type //is_void<U> is always false for us, used to make the deduction argument and SFINAE work
                evaluate_f(U lambda)
                {
                        return g(lambda);
                }

                //! Function version of a std::function () operator, use SFINAE to differentiate from the two possibilities of instantiation, not derived case
                /*!
                 \tparam U used for SFINAE purpose, is the type of lambda parameter, should not be void
                 \param lambda the independent variable for function evaluation, usually a Real or a VectorXr
                 \return a Ctype: evaluation of the function at point lambda
                 \note never choose U void (it is unresonable and beaks the mechanism)
                */
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)!=0 || std::is_void<U>::value, Ctype>::type
                evaluate_f(U lambda )
                {
                        return this->compute_f(lambda);
                }

                //! Evaluation of first derivative, if derived
                /*!
                 \tparam U used for SFINAE purpose, is the type of lambda parameter, should not be void
                 \param lambda the independent variable for function evaluation, usually a Real or a VectorXr
                 \return a Tuple: evaluation of the first derivative at point lambda
                 \note never choose U void (it is unresonable and beaks the mechanism)
                */
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)==0 || std::is_void<U>::value, Tuple>::type //is_void<U> is always false for us, used to make the deduction argument and SFINAE work
                evaluate_first_derivative(U lambda)
                {
                        return dg(lambda);
                }

                //! Evaluation of first derivative, if not derived
                /*!
                 \tparam U used for SFINAE purpose, is the type of lambda parameter, should not be void
                 \param lambda the independent variable for function evaluation, usually a Real or a VectorXr
                 \return a Tuple: evaluation of the first derivative at point lambda
                 \note never choose U void (it is unresonable and beaks the mechanism)
                */
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)!=0 || std::is_void<U>::value, Tuple>::type
                evaluate_first_derivative(U lambda)
                {
                        return this->compute_fp(lambda);
                }

                //! Evaluation of second derivative, if derived
                /*!
                 \tparam U used for SFINAE purpose, is the type of lambda parameter, should not be void
                 \param lambda the independent variable for function evaluation, usually a Real or a VectorXr
                 \return a Hessian: evaluation of the second derivative at point lambda
                 \note never choose U void (it is unresonable and beaks the mechanism)
                */
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)==0 || std::is_void<U>::value, Hessian>::type //is_void<U> is always false for us, used to make the deduction argument and SFINAE work
                evaluate_second_derivative(U lambda)
                {
                        return ddg(lambda);
                }

                //! Evaluation of second derivative, if not derived
                /*!
                 \tparam U used for SFINAE purpose, is the type of lambda parameter, should not be void
                 \param lambda the independent variable for function evaluation, usually a Real or a VectorXr
                 \return a Hessian: evaluation of the second derivative at point lambda
                 \note never choose U void (it is unresonable and beaks the mechanism)
                */
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)!=0 || std::is_void<U>::value, Hessian>::type
                evaluate_second_derivative(U lambda )
                {
                        return this->compute_fs(lambda);
                }
};

#endif
