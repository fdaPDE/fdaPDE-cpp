#ifndef _FPIRLS_H
#define _FPIRLS_H

#include <cmath>
#include <math.h>
#include <array>

#include "mixedFERegression.h"
#include "evaluator.h"
#include "fdaPDE.h"

#include"printer.h"


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Base {

  protected:

  const MeshHandler<ORDER, mydim, ndim> &mesh_;
  InputHandler inputData_; // it contains the original data of the problem
  MixedFERegression<InputHandler, Integrator,ORDER, IntegratorGaussP3, 0, 0, mydim, ndim>  regression_;

  //RegressionData& pseudoData_; // contains data used to perform step (2) of f-PIRLS algo
  // question is: should I introduce a new template for this (data in step (2))? At the moment it is RegressionData, but one can imagine future developements around childs of regressionData...

  std::vector<VectorXr> mu_; // mean vector
  std::vector<VectorXr> pseudoObservations_; // pseudo observations
  std::vector<VectorXr> G_; // diag(link_deriv(mu)) it is a vector since it would be more memory consuming to declere it as a matrix
  std::vector<VectorXr> WeightsMatrix_; // (G^-2 * Var_function(mu)^-1) it is a vector for the same reason as above

  std::vector<std::array<Real,2>> current_J_values;
  std::vector<std::array<Real,2>> past_J_values; // stores the value of the functional J at each iteration in order to apply the stopping criterion
  // the value of the functional is saved deparated (parametric and non-parametric part)

  std::vector<UInt> n_iterations; // current n° of iteration of PIRLS

  VectorXr forcingTerm;
  bool isSpaceVarying = false;

  MatrixXv _solution; //!A Eigen::VectorXr: Stores the system solution.
  MatrixXr _dof; //! A vector storing the computed dofs
  
  std::vector<Real> _GCV; //! A vector storing GCV values
  std::vector<Real> _J_minima;

  // Evaluation of the solution in the locations and beta estimates
  MatrixXv _beta_hat;
  MatrixXv _fn_hat;

  UInt bestLambdaS_ = 0;  //!Stores the index of the best lambdaS according to GCV
  Real _bestGCV = 10e20;  //!Stores the value of the best GCV




    void compute_pseudoObs(UInt& lambda_index); // perform step (1) of f-PIRLS

    void compute_G(UInt& lambda_index); //assemble G matrix

    void compute_Weights(UInt& lambda_index); // assemble W matrix

    void update_solution(UInt& lambda_index); // perform step (2) of f-PIRLS it is dependent on templates specs

    void compute_mu(UInt& lambda_index); // perform step (3) of f-PIRLS

    bool stopping_criterion(UInt& lambda_index); // it stops PIRLS based on difference between functionals J_k J_k+1 or n_iterations > max_num_iterations

    std::array<Real,2> compute_J(UInt& lambda_index); // compute the current value of J

    void compute_GCV(UInt& lambda_index); // compute the GCV value for a given lambda

    // link and other functions, it can be extended to handle vector of parameters
    virtual Real link(const Real& mu)const = 0; // g(.)

    virtual Real link_deriv(const Real& mu)const = 0; // g'(.)

    virtual Real inv_link(const Real& theta)const = 0; // g^-1(.)

    virtual Real var_function(const Real& mu)const = 0; // V(mu)

    virtual Real dev_function(const Real&mu, const Real& x)const = 0; //deviation function: used as norm in GCV


  public:

    FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0); // Constructor

    //! A destructor
   virtual ~FPIRLS_Base(){};

   // Main method: perform PIRLS and instanciate the solution in _solution , _dof
   void apply(const ForcingTerm& u);

   //! A inline member that returns a VectorXr, returns the whole solution_.
   inline MatrixXv const & getSolution() const{return _solution;}
   //! A function returning the computed dofs of the model
   inline MatrixXr const & getDOF() const{return _dof;}
   //! A function returning the current value of J
   inline std::vector<Real> const & get_J() const{return _J_minima;}
   //! A inline member that returns a VectorXr, returns the final beta estimate.
   inline MatrixXv const & getBetaEst() const{return _beta_hat;}
   //! A inline member that returns a VectorXr, returns the final function estimates.
   inline MatrixXv const & getFunctionEst() const{return _fn_hat;}
   //! A inline member that returns scale parameter estimates.
   inline std::vector<Real> const getScaleParamEst() const{return std::vector<Real>(1,1) ;}
   //! A inline member that returns a the computed (or not) GCV estimates. If GCV is not computed, -1 is returned
   inline std::vector<Real> const & getGCV() const{return _GCV;}
   //! A method returning the index of the best lambdaS according to GCV
   inline UInt getBestLambdaS(){return bestLambdaS_;}
   

    //! A function returning the computed barycenters of the locationss
   inline MatrixXr const & getBarycenters() const{return regression_.getBarycenters();}
   //! A function returning the element ids of the locations
   inline VectorXi const & getElementIds() const{return regression_.getElementIds();}
   
   // FURTHER QUESTIONS - where to implement GCV (lambda can be specified by user or estimated via GCV), where to implement the estimate of phi

   // QUESTION: se avessi un costruttore padre con parametri e non definisco il costruttore figlio (aka, figlio dispone solo di default constr.) e creo un oggetto figlio chiamando il costruttore con parametri ottengo automaticamente quello del padre? In caso negativo devo implementare il costruttore anche per le singole distribuzioni.

};

//Template Specialization
// Laplace or Elliptic
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS: public FPIRLS_Base< InputHandler,  Integrator,  ORDER,  mydim,  ndim>{

  public:

    FPIRLS(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS_Base<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0){};
     //! A destructor
   virtual ~FPIRLS(){};

   virtual void apply();
};

//EllipticSpaceVarying
template <typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS< GAMDataEllipticSpaceVarying,  Integrator,  ORDER,  mydim,  ndim>: public FPIRLS_Base< GAMDataEllipticSpaceVarying,  Integrator,  ORDER,  mydim,  ndim>{

  public:

    FPIRLS(const MeshHandler<ORDER,mydim,ndim>& mesh, GAMDataEllipticSpaceVarying& inputData, VectorXr mu0):
      FPIRLS_Base<GAMDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0){};

    virtual ~FPIRLS(){};

    virtual void apply();
};



// classes of distribution having scale params
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Scaled: public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

  protected:

    Real scale_parameter_;
    bool scale_parameter_flag_;
    std::vector<Real> _scale_parameter_estimates;

    void compute_scale_param(); //compute the estimates of the scale param

  public:

  FPIRLS_Scaled(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, Real scale_parameter, bool scale_parameter_flag):
    FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0), scale_parameter_(scale_parameter), scale_parameter_flag_(scale_parameter_flag){};

  //! A destructor
  virtual ~FPIRLS_Scaled(){};

  //! Overriden apply method taking into account also the estimates of scale parameter
  void apply() override;

  //! A inline member that returns scale parameter estimates.
  inline std::vector<Real> const getScaleParamEst() const{return _scale_parameter_estimates;}

};


//------------- Family Distributions Spcecification ----------------

// Bernoulli
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Bernoulli : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

  protected:
    inline Real link(const Real& mu)const{ return log(mu/(1 - mu)); }

    inline Real inv_link(const Real& theta)const{ return 1/(1 + exp(-theta)); }

    inline Real link_deriv(const Real& mu)const{ return 1/(mu*(1-mu)); }

    inline Real var_function(const Real& mu)const{ return(mu*(1-mu)); }

    inline Real dev_function(const Real& mu, const Real& x)const{return (x == 0)? 2*log(1/(1-mu)) : 2*log(1/mu);}

  //  void initialize_mu(const VectorXr& y);

  public:

    FPIRLS_Bernoulli(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0){};
};


// Probit
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Probit : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {
  private: 
    Real erf_inv(const Real& x)const;
  
  protected:
    inline Real link(const Real& mu)const{ return erf_inv(mu); }

    inline Real inv_link(const Real& theta)const{ return erf(theta); }

    inline Real link_deriv(const Real& mu)const{ return 1/(sqrt(2*M_PI)) * exp(- 0.5*(mu*mu)) ;}

    inline Real var_function(const Real& mu)const{ return (mu*(1-mu)); }

    inline Real dev_function(const Real& mu, const Real& x)const{return (x == 0)? 2*log(1/(1-mu)) : 2*log(1/mu);}

  //  void initialize_mu(const VectorXr& y);

  public:

    FPIRLS_Probit(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0){};
};


// cLogLog
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_cLogLog : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

  protected:
    inline Real link(const Real& mu)const{ return log(-log(1-mu)); }

    inline Real inv_link(const Real& theta)const{ return  1-exp(-exp(theta)) ; }

    inline Real link_deriv(const Real& mu)const{ return  -1/(log(1-mu))*1/(1-mu); }

    inline Real var_function(const Real& mu)const{ return(mu*(1-mu)); }

    inline Real dev_function(const Real& mu, const Real& x)const{return (x == 0)? 2*log(1/(1-mu)) : 2*log(1/mu);}

  //  void initialize_mu(const VectorXr& y);

  public:

    FPIRLS_cLogLog(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0){};
};




// Poisson
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
    class FPIRLS_Poisson : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu)const{ return log(mu); }

      inline Real link_deriv(const Real& mu)const{ return 1/mu; }

      inline Real inv_link(const Real& theta)const{ return exp(theta); }

      inline Real var_function(const Real& mu)const{ return mu ;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return (x>0) ? x*log(x/mu) - (x-mu): mu; }
    /*  void initialize_mu(const VectorXr & y) {
      this->mu_ = y;} // It is different for binary or non-binary outcomes
    */
    public:

    FPIRLS_Poisson(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0){};

};

// Exponential
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
    class FPIRLS_Exponential : public FPIRLS <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu)const{ return -1/mu; }

      inline Real link_deriv(const Real& mu)const{ return 1/(mu*mu); }

      inline Real inv_link(const Real& theta)const{ return - 1/theta; }

      inline Real var_function(const Real& mu)const{ return mu*mu ;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return 2*(((x-mu)/mu)-log(x/mu)); }

    public:

    FPIRLS_Exponential(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):
      FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0){};

};

//------------- Scaled Distributions ----------

// Gamma
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
    class FPIRLS_Gamma : public FPIRLS_Scaled <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu)const{ return - 1/mu ; }

      inline Real link_deriv(const Real& mu)const{ return 1/(mu*mu); }

      inline Real inv_link(const Real& theta)const{ return - 1/theta; }

      inline Real var_function(const Real& mu)const{ return mu*mu ;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return 2*(((x-mu)/mu)-log(x/mu)); }

    /*  void initialize_mu(const VectorXr & y) {
      this->mu_ = y;} // It is different for binary or non-binary outcomes
    */
    public:

    FPIRLS_Gamma(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, Real scale_parameter,bool scale_parameter_flag):
      FPIRLS_Scaled<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, scale_parameter, scale_parameter_flag){};

};

// Gaussian
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
    class FPIRLS_Gaussian : public FPIRLS_Scaled <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu)const{ return mu ; }

      inline Real link_deriv(const Real& mu)const{ return 1; }

      inline Real inv_link(const Real& theta)const{ return theta; }

      inline Real var_function(const Real& mu)const{ return 1 ;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return (x-mu)*(x-mu);}

    /*  void initialize_mu(const VectorXr & y) {
      this->mu_ = y;} // It is different for binary or non-binary outcomes
    */
    public:

    FPIRLS_Gaussian(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, Real scale_parameter,bool scale_parameter_flag):
      FPIRLS_Scaled<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, scale_parameter, scale_parameter_flag){};

};


// InvGaussian
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
    class FPIRLS_InvGaussian : public FPIRLS_Scaled <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu)const{ return -1/(2*mu*mu) ; }

      inline Real link_deriv(const Real& mu)const{ return 1/(mu*mu*mu)   ; }

      inline Real inv_link(const Real& theta)const{ return sqrt(-1/(2*theta)); }

      inline Real var_function(const Real& mu)const{ return mu*mu*mu;}

      inline Real dev_function(const Real&mu, const Real& x)const{ return (x-mu)*(x-mu)/(mu*mu*x);}

    /*  void initialize_mu(const VectorXr & y) {
      this->mu_ = y;} // It is different for binary or non-binary outcomes
    */
    public:

    FPIRLS_InvGaussian(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0, Real scale_parameter,bool scale_parameter_flag):
      FPIRLS_Scaled<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData, mu0, scale_parameter, scale_parameter_flag){};

};

// Probit
// TO BE DONE: how to get normal cdf and its inverse? thay're needed for probit

// Cloglog
// TO BE DONE: W. write his wrong in the old R code, has to be fixed here




#include "FPIRLS_imp.h"

#endif
