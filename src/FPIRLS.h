#ifndef _FPIRLS_H
#define _FPIRLS_H

#include <cmath>
#include <math.h>
#include <array>

#include "mixedFERegression.h"
#include "evaluator.h"
#include "fdaPDE.h"


//! @brief An abstract class that implement the apply method for the FPIRLS algorithm.
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Base {

  protected:

   const MeshHandler<ORDER, mydim, ndim> &mesh_;
   InputHandler& inputData_; //!< It contains the data of the problem (RegressionDataGAM)
   MixedFERegression<InputHandler, Integrator,ORDER, IntegratorGaussP3, 0, 0, mydim, ndim>  regression_;

  
   std::vector<VectorXr> mu_; //!< Mean vector
   std::vector<VectorXr> pseudoObservations_; //! Pseudodata observations
   std::vector<VectorXr> G_; //!< diag(link_deriv(mu)) it is a vector since it would be more memory consuming to declere it as a matrix
   std::vector<VectorXr> WeightsMatrix_; //!< (G^-2 * Var_function(mu)^-1) it is a vector because it is a diagonal matrix

   // the value of the functional is saved deparated (parametric and non-parametric part)
   std::vector<std::array<Real,2>> current_J_values;
   std::vector<std::array<Real,2>> past_J_values; //!< Stores the value of the functional J at each iteration in order to apply the stopping criterion

   std::vector<UInt> n_iterations; //!< Current number of iteration of PIRLS

   VectorXr forcingTerm;
   bool isSpaceVarying = false; //!< True only in space varying case.

   MatrixXv _solution; //!< Stores the system solution.
   MatrixXr _dof; //!< A matrix of VectorXr storing the computed dofs.
  
   std::vector<Real> _GCV; //!< A vector storing GCV values.
   std::vector<Real> _J_minima;//!< Stores the minimum value of the functional.

   // Evaluation of the solution in the locations and beta estimates
   MatrixXv _beta_hat;//!< Betas estimated if the model has covariates.
   MatrixXv _fn_hat; //!< Function coefficients estimated.

   UInt bestLambdaS_ = 0;  //!< Stores the index of the best lambdaS according to GCV.
   Real _bestGCV = 10e20;  //!< Stores the value of the best GCV.

   bool scale_parameter_flag_; //!< True if the distribution has the scale parameter and if it is not a given input.
   Real _scale_param;
   std::vector<Real> _variance_estimates; //!< Stores the variance estimates for each lambda.

   //! A method that computes the pseudodata. It perform step (1) of f-PIRLS.
   void compute_pseudoObs(UInt& lambda_index);
   //! A method that assembles G matrix.
   void compute_G(UInt& lambda_index);
   //! A method that assembles the weights matrix ( a diagonal matrix, hence it is stored as vector).
   void compute_Weights(UInt& lambda_index);
   //! A method that updates the solution. It perform step (2) of F-PIRLS.
   void update_solution(UInt& lambda_index);
   //! A method that updates mu vector. It perform step (3) of F-PIRLS.
   void compute_mu(UInt& lambda_index);
   //! A method that stops PIRLS based on difference between functionals J_k J_k+1 or n_iterations > max_num_iterations .
   bool stopping_criterion(UInt& lambda_index);
   //! A method that computes and return the current value of the functional J. It is divided in parametric and non parametric part.
   std::array<Real,2> compute_J(UInt& lambda_index);
   //! A method that computes the GCV value for a given lambda.
   void compute_GCV(UInt& lambda_index);
   //! A method that computes the estimates of the variance. It depends on the scale flags: only the Gamma and InvGaussian distributions have the scale parameter.
   void compute_variance_est();

   // link and other functions. Definited as pure virtual methods, the implementaton depend on the choosen distributions 

   //! A pure virtual method that represents the link function: g(.)
   virtual Real link(const Real& mu)const = 0; 
   //! A pure virtual method that represents the derivative function: g'(.)
   virtual Real link_deriv(const Real& mu)const = 0;
   //! A pure virtual method that represents the inverse link function: g^-1(.)
   virtual Real inv_link(const Real& theta)const = 0;
   //! A pure virtual method that represents the variance function: V(mu)
   virtual Real var_function(const Real& mu)const = 0;
   //! A pure virtual method that represents the deviation function: used as norm in GCV
   virtual Real dev_function(const Real&mu, const Real& x)const = 0; 


  public:

    FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, bool scale_parameter_flag, Real scale_param); // Constructor

    //! A virutal destructor
   virtual ~FPIRLS_Base(){};

   //! Main method: perform PIRLS and instanciate the solution in _solution , _dof
   void apply(const ForcingTerm& u);

   //! An inline member that returns a VectorXr, returns the whole solution_.
   inline MatrixXv const & getSolution() const{return _solution;}
   //! A method returning the computed dofs of the model
   inline MatrixXr const & getDOF() const{return _dof;}
   //! A method returning the current value of J
   inline std::vector<Real> const & get_J() const{return _J_minima;}
   //! A inline member that returns a VectorXr, returns the final beta estimate.
   inline MatrixXv const & getBetaEst() const{return _beta_hat;}
   //! An inline member that returns a VectorXr, returns the final function estimates.
   inline MatrixXv const & getFunctionEst() const{return _fn_hat;}
   //! An inline member that returns the variance estimates.
   inline std::vector<Real> const getVarianceEst() const{return  _variance_estimates;}
   //! An inline member that returns a the computed (or not) GCV estimates. If GCV is not computed, -1 is returned
   inline std::vector<Real> const & getGCV() const{return _GCV;}
   //! A method returning the index of the best lambdaS according to GCV
   inline UInt getBestLambdaS(){return bestLambdaS_;}

   
   //! A method returning the computed barycenters of the locationss
   inline MatrixXr const & getBarycenters() const{return regression_.getBarycenters();}
   //! A method returning the element ids of the locations
   inline VectorXi const & getElementIds() const{return regression_.getElementIds();}
   

};


//------------- Template Specialization ----------------


//! @brief A class used for the FPIRLS_base template specialization: Laplace or Elliptic cases.
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS: public FPIRLS_Base< InputHandler,  Integrator,  ORDER,  mydim,  ndim>{

  public:

    FPIRLS(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, bool scale_parameter_flag, Real scale_param):
      FPIRLS_Base<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, scale_parameter_flag, scale_param){};
     //! A virtual destructor
   virtual ~FPIRLS(){};

   virtual void apply();
};

//! @brief A class used for the FPIRLS_base template specialization: EllipticSpaceVarying case.
template <typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS< GAMDataEllipticSpaceVarying,  Integrator,  ORDER,  mydim,  ndim>: public FPIRLS_Base< GAMDataEllipticSpaceVarying,  Integrator,  ORDER,  mydim,  ndim>{

  public:

    FPIRLS(const MeshHandler<ORDER,mydim,ndim>& mesh, GAMDataEllipticSpaceVarying& inputData, VectorXr mu0, bool scale_parameter_flag, Real scale_param):
      FPIRLS_Base<GAMDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, scale_parameter_flag, scale_param){};

    virtual ~FPIRLS(){};

    virtual void apply();
};



//------------- Family Distributions Spcecification ----------------

//! @brief A class that specify the Bernoulli distribution for the FPIRLS class.
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Bernoulli : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

  protected:
    inline Real link(const Real& mu)const{ return log(mu/(1 - mu)); }

    inline Real inv_link(const Real& theta)const{ return 1/(1 + exp(-theta)); }

    inline Real link_deriv(const Real& mu)const{ return 1/(mu*(1-mu)); }

    inline Real var_function(const Real& mu)const{ return(mu*(1-mu)); }

    inline Real dev_function(const Real& mu, const Real& x)const{return (x == 0)? 2*log(1/(1-mu)) : 2*log(1/mu);}

  public:

    FPIRLS_Bernoulli(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, false, 1){};
};


//! @brief A class that specify the Poisson distribution for the FPIRLS class.
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Poisson : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu)const{ return log(mu); }

      inline Real link_deriv(const Real& mu)const{ return 1/mu; }

      inline Real inv_link(const Real& theta)const{ return exp(theta); }

      inline Real var_function(const Real& mu)const{ return mu ;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return (x>0) ? x*log(x/mu) - (x-mu): mu; }

    public:

    FPIRLS_Poisson(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, false, 1){};

};

//! @brief A class that specify the Exponential distribution for the FPIRLS class.
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Exponential : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> 
{

    protected:

      inline Real link(const Real& mu)const{ return -1/mu; }

      inline Real link_deriv(const Real& mu)const{ return 1/(mu*mu); }

      inline Real inv_link(const Real& theta)const{ return - 1/theta; }

      inline Real var_function(const Real& mu)const{ return mu*mu ;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return 2*(((x-mu)/mu)-log(x/mu)); }

    public:

    FPIRLS_Exponential(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, false, 1){};

};

//------------- Scaled Distributions ----------

//! @brief A class that specify the Gamma distribution for the FPIRLS class.
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Gamma : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu)const{ return - 1/mu ; }

      inline Real link_deriv(const Real& mu)const{ return 1/(mu*mu); }

      inline Real inv_link(const Real& theta)const{ return - 1/theta; }

      inline Real var_function(const Real& mu)const{ return mu*mu ;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return 2*(((x-mu)/mu)-log(x/mu)); }

    public:

    FPIRLS_Gamma(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, bool scale_parameter_flag, Real scale_param):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, scale_parameter_flag, scale_param){};

};



#include "FPIRLS_imp.h"

#endif
