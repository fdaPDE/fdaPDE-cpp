// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/core/utils/Symbols.h>
#include <fdaPDE/models/regression/SRPDE.h>
using fdaPDE::models::SRPDE;
#include <fdaPDE/core/utils/DataStructures/BlockFrame.h>
#include <fdaPDE/models/ModelTraits.h>
using fdaPDE::models::Sampling;
#include <fdaPDE/core/FEM/PDE.h>
using fdaPDE::core::FEM::DefaultOperator;
using fdaPDE::core::FEM::PDE;
#include <fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h>
using fdaPDE::core::FEM::SpaceVaryingAdvection;
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
using fdaPDE::core::FEM::SpaceVaryingReaction;
#include <fdaPDE/core/MESH/Mesh.h>
using fdaPDE::core::MESH::Mesh;
#include <fdaPDE/models/SamplingDesign.h>
using fdaPDE::models::Sampling;

// this file contains the R wrapper for the SRPDE model

template <unsigned int M, unsigned int N, unsigned int R, typename F>
class RegularizingPDE {
private:
  typedef typename std::decay<F>::type BilinearFormType;
  // internal data
  Mesh<M,N,R> domain_;
  PDE<M,N,R, BilinearFormType, DMatrix<double>> pde_;
public:
  // constructor
  RegularizingPDE(const Rcpp::List& R_Mesh) :
    // initialize domain
    domain_(Rcpp::as<DMatrix<double>>(R_Mesh["nodes"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["edges"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["elements"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["neigh"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["boundary"])),
    pde_(domain_) { pde_.setBilinearForm(BilinearFormType()); };
  
  // setters
  void set_dirichlet_bc(const DMatrix<double>& data){ pde_.setDirichletBC(data); }
  void set_forcing_term(const DMatrix<double>& data){ pde_.setForcing(data); }
  // getters
  DMatrix<double> get_quadrature_nodes() const { return pde_.integrator().quadratureNodes(domain_); };
  DMatrix<double> get_dofs_coordinates() const { return domain_.dofCoords(); };
  const PDE<M,N,R, BilinearFormType, DMatrix<double>>& pde() { return pde_; }
  
  // compile time informations
  typedef PDE<M,N,R, BilinearFormType, DMatrix<double>> PDEType;
};
// define 2D simple Laplacian regularization.
typedef RegularizingPDE<2,2,1, decltype( std::declval<Laplacian<DefaultOperator>>() )>
Laplacian_2D_Order1;
// expose RegularizingPDE as possible argument to other Rcpp modules
RCPP_EXPOSED_AS  (Laplacian_2D_Order1)
RCPP_EXPOSED_WRAP(Laplacian_2D_Order1)

RCPP_MODULE(Laplacian_2D_Order1) {
  Rcpp::class_<Laplacian_2D_Order1>("Laplacian_2D_Order1")
    .constructor<Rcpp::List>()
    // getters
    .method("get_quadrature_nodes", &Laplacian_2D_Order1::get_quadrature_nodes)
    .method("get_dofs_coordinates", &Laplacian_2D_Order1::get_dofs_coordinates)
    // setters
    .method("set_dirichlet_bc",     &Laplacian_2D_Order1::set_dirichlet_bc)
    .method("set_forcing_term",     &Laplacian_2D_Order1::set_forcing_term);
}

// constant coefficients PDE type
template<unsigned int M>
using ConstantCoefficientsPDE =
  decltype( std::declval<Laplacian<SMatrix<M>>>() +
	    std::declval<Gradient <SVector<M>>>() +
	    std::declval<Identity <double>>() );

// specialization of RegularizingPDE for constant coefficients case
template <unsigned int M, unsigned int N, unsigned int R>
class RegularizingPDE<M,N,R, ConstantCoefficientsPDE<M>> {
private:
  typedef typename std::decay<ConstantCoefficientsPDE<M>>::type BilinearFormType;
  // internal data
  Mesh<M,N,R> domain_;
  PDE<M,N,R, BilinearFormType, DMatrix<double>> pde_;
public:
  // constructor
  RegularizingPDE(const Rcpp::List& R_Mesh) :
    // initialize domain
    domain_(Rcpp::as<DMatrix<double>>(R_Mesh["nodes"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["edges"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["elements"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["neigh"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["boundary"])),
    pde_(domain_) {};
  
  // setters
  void set_dirichlet_bc(const DMatrix<double>& data){ pde_.setDirichletBC(data); }
  void set_forcing_term(const DMatrix<double>& data){ pde_.setForcing(data); }
  void set_PDE_parameters(const Rcpp::List& data){
    SMatrix<M> K = Rcpp::as<DMatrix<double>>(data["diffusion"]);
    SVector<M> b = Rcpp::as<DVector<double>>(data["transport"]);
    double c = Rcpp::as<double>(data["reaction"]);
    BilinearFormType bilinearForm = Laplacian(K) + Gradient(b) + Identity(c);
    pde_.setBilinearForm(bilinearForm);  
  };

  // getters
  DMatrix<double> get_quadrature_nodes() const { return pde_.integrator().quadratureNodes(domain_); };
  DMatrix<double> get_dofs_coordinates() const { return domain_.dofCoords(); };
  const PDE<M,N,R, BilinearFormType, DMatrix<double>>& pde() { return pde_; }
  
  // compile time informations
  typedef PDE<M,N,R, BilinearFormType, DMatrix<double>> PDEType;
};

// define 2D costant coefficient PDE regularization.
typedef RegularizingPDE<2,2,1, ConstantCoefficientsPDE<2>> ConstantCoefficients_2D_Order1;
// expose RegularizingPDE as possible argument to other Rcpp modules
RCPP_EXPOSED_AS  (ConstantCoefficients_2D_Order1)
RCPP_EXPOSED_WRAP(ConstantCoefficients_2D_Order1)

RCPP_MODULE(ConstantCoefficients_2D_Order1) {
  Rcpp::class_<ConstantCoefficients_2D_Order1>("ConstantCoefficients_2D_Order1")
    .constructor<Rcpp::List>()
    // getters
    .method("get_quadrature_nodes", &ConstantCoefficients_2D_Order1::get_quadrature_nodes)
    .method("get_dofs_coordinates", &ConstantCoefficients_2D_Order1::get_dofs_coordinates)
    // setters
    .method("set_dirichlet_bc",     &ConstantCoefficients_2D_Order1::set_dirichlet_bc)
    .method("set_forcing_term",     &ConstantCoefficients_2D_Order1::set_forcing_term)
    .method("set_PDE_parameters",   &ConstantCoefficients_2D_Order1::set_PDE_parameters);
}

// space varying PDE type
template<unsigned int M>
using SpaceVaryingPDE =
  decltype( std::declval< Laplacian< decltype(std::declval<SpaceVaryingDiffusion<M>>().asParameter()) >>() +
	    std::declval< Gradient < decltype(std::declval<SpaceVaryingAdvection<M>>().asParameter()) >>() +
	    std::declval< Identity < decltype(std::declval<SpaceVaryingReaction>().asParameter()) >>() );

// specialization of RegularizingPDE for space varying
template <unsigned int M, unsigned int N, unsigned int R>
class RegularizingPDE<M,N,R, SpaceVaryingPDE<M>> {
private:
  typedef typename std::decay<SpaceVaryingPDE<M>>::type BilinearFormType;
  // internal data
  Mesh<M,N,R> domain_;
  PDE<M,N,R, BilinearFormType, DMatrix<double>> pde_;
  // space-varying functors
  SpaceVaryingDiffusion<M> diffusion_;
  SpaceVaryingAdvection<M> advection_;
  SpaceVaryingReaction     reaction_;
public:
  // constructor
  RegularizingPDE(const Rcpp::List& R_Mesh) :
    // initialize domain
    domain_(Rcpp::as<DMatrix<double>>(R_Mesh["nodes"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["edges"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["elements"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["neigh"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["boundary"])),
    pde_(domain_) {};
  
  // setters
  void set_dirichlet_bc(const DMatrix<double>& data){ pde_.setDirichletBC(data); }
  void set_forcing_term(const DMatrix<double>& data){ pde_.setForcing(data); }
  void set_PDE_parameters(const Rcpp::List& data){
    DMatrix<double> K = Rcpp::as<DMatrix<double>>(data["diffusion"]);
    diffusion_.setData(K);
    DMatrix<double> b = Rcpp::as<DMatrix<double>>(data["transport"]);
    advection_.setData(b);
    DMatrix<double> c = Rcpp::as<DMatrix<double>>(data["reaction"]);
    reaction_.setData(c);
    BilinearFormType bilinearForm = Laplacian(diffusion_.asParameter()) + Gradient(advection_.asParameter()) + Identity(reaction_.asParameter());
    pde_.setBilinearForm(bilinearForm);  
  };

  // getters
  DMatrix<double> get_quadrature_nodes() const { return pde_.integrator().quadratureNodes(domain_); };
  DMatrix<double> get_dofs_coordinates() const { return domain_.dofCoords(); };
  const PDE<M,N,R, BilinearFormType, DMatrix<double>>& pde() { return pde_; }
  
  // compile time informations
  typedef PDE<M,N,R, BilinearFormType, DMatrix<double>> PDEType;
};

// define 2D costant coefficient PDE regularization.
typedef RegularizingPDE<2,2,1, SpaceVaryingPDE<2>> SpaceVarying_2D_Order1;
// expose RegularizingPDE as possible argument to other Rcpp modules
RCPP_EXPOSED_AS  (SpaceVarying_2D_Order1)
RCPP_EXPOSED_WRAP(SpaceVarying_2D_Order1)

RCPP_MODULE(SpaceVarying_2D_Order1) {
  Rcpp::class_<SpaceVarying_2D_Order1>("SpaceVarying_2D_Order1")
    .constructor<Rcpp::List>()
    // getters
    .method("get_quadrature_nodes", &SpaceVarying_2D_Order1::get_quadrature_nodes)
    .method("get_dofs_coordinates", &SpaceVarying_2D_Order1::get_dofs_coordinates)
    // setters
    .method("set_dirichlet_bc",     &SpaceVarying_2D_Order1::set_dirichlet_bc)
    .method("set_forcing_term",     &SpaceVarying_2D_Order1::set_forcing_term)
    .method("set_PDE_parameters",   &SpaceVarying_2D_Order1::set_PDE_parameters);
}

// wrapper for SRPDE module
template <typename RegularizingPDE_, Sampling S_> class R_SRPDE;

// macro for the definition of common functionalities to all SRPDE Rcpp modules
#define SRPDE_MODULE(REGULARIZATION_TYPE, SAMPLING_TYPE, ... )                 \
  template <> class R_SRPDE<REGULARIZATION_TYPE, SAMPLING_TYPE> {              \
  protected:                                                                   \
    typedef REGULARIZATION_TYPE RegularizingPDE_;                              \
    RegularizingPDE_ regularization_;                                          \
    /* the model this Rcpp module wraps */                                     \
    SRPDE<typename RegularizingPDE_::PDEType, SAMPLING_TYPE> model_;           \
    BlockFrame<double, int> df_;                                               \
                                                                               \
  public:                                                                      \
    /* constructor */                                                          \
    R_SRPDE(const RegularizingPDE_ &regularization)                            \
        : regularization_(regularization) {                                    \
      model_.setPDE(regularization_.pde());                                    \
    };                                                                         \
                                                                               \
    /* setters */                                                              \
    void set_lambda_s(double lambdaS) { model_.setLambdaS(lambdaS); }          \
    void set_observations(const DMatrix<double> &data) {                       \
      df_.template insert<double>(OBSERVATIONS_BLK, data);                     \
    }                                                                          \
    void set_covariates(const DMatrix<double> &data) {                         \
      df_.template insert<double>(DESIGN_MATRIX_BLK, data);                    \
    }                                                                          \
    /* getters */                                                              \
    DMatrix<double> f() const { return model_.f(); }                           \
    DMatrix<double> beta() const { return model_.beta(); }                     \
                                                                               \
    /* initialize model and solve smoothing problem */                         \
    void solve() {                                                             \
      model_.setData(df_);                                                     \
      model_.init();                                                           \
      model_.solve();                                                          \
    }                                                                          \
    /* additional methods */						       \
    __VA_ARGS__								       \
  };

// Laplacian regularization, sampling at mesh nodes
SRPDE_MODULE(Laplacian_2D_Order1, Sampling::GeoStatMeshNodes); 
// definition of Rcpp module
typedef R_SRPDE<Laplacian_2D_Order1, Sampling::GeoStatMeshNodes>
SRPDE_Laplacian_2D_GeoStatNodes;
RCPP_MODULE(SRPDE_Laplacian_2D_GeoStatNodes) {
  Rcpp::class_<SRPDE_Laplacian_2D_GeoStatNodes>("SRPDE_Laplacian_2D_GeoStatNodes")
    .constructor<Laplacian_2D_Order1>()
    // getters
    .method("f",                &SRPDE_Laplacian_2D_GeoStatNodes::f)
    .method("beta",             &SRPDE_Laplacian_2D_GeoStatNodes::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_Laplacian_2D_GeoStatNodes::set_lambda_s)
    .method("set_observations", &SRPDE_Laplacian_2D_GeoStatNodes::set_observations)
    .method("set_covariates",   &SRPDE_Laplacian_2D_GeoStatNodes::set_covariates)
    // solve method
    .method("solve",            &SRPDE_Laplacian_2D_GeoStatNodes::solve);
}

// Laplacian regularization, sampling at given locations
SRPDE_MODULE(Laplacian_2D_Order1, Sampling::GeoStatLocations,
	     void set_locations(const DMatrix<double>& data) { model_.setLocs(data); }
	     )
// definition of Rcpp module
typedef R_SRPDE<Laplacian_2D_Order1, Sampling::GeoStatLocations>
SRPDE_Laplacian_2D_GeoStatLocations;
RCPP_MODULE(SRPDE_Laplacian_2D_GeoStatLocations) {
  Rcpp::class_<SRPDE_Laplacian_2D_GeoStatLocations>("SRPDE_Laplacian_2D_GeoStatLocations")
    .constructor<Laplacian_2D_Order1>()
    // getters
    .method("f",                &SRPDE_Laplacian_2D_GeoStatLocations::f)
    .method("beta",             &SRPDE_Laplacian_2D_GeoStatLocations::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_Laplacian_2D_GeoStatLocations::set_lambda_s)
    .method("set_observations", &SRPDE_Laplacian_2D_GeoStatLocations::set_observations)
    .method("set_covariates",   &SRPDE_Laplacian_2D_GeoStatLocations::set_covariates)
    .method("set_locations",    &SRPDE_Laplacian_2D_GeoStatLocations::set_locations)
    // solve method
    .method("solve",            &SRPDE_Laplacian_2D_GeoStatLocations::solve);
}

// Laplacian regularization, areal sampling
SRPDE_MODULE(Laplacian_2D_Order1, Sampling::Areal,
	     void set_subdomains(const DMatrix<int>& data) { model_.setSubdomains(data); }
	     )
// definition of Rcpp module
typedef R_SRPDE<Laplacian_2D_Order1, Sampling::Areal>
SRPDE_Laplacian_2D_Areal;
RCPP_MODULE(SRPDE_Laplacian_2D_Areal) {
  Rcpp::class_<SRPDE_Laplacian_2D_Areal>("SRPDE_Laplacian_2D_Areal")
    .constructor<Laplacian_2D_Order1>()
    // getters
    .method("f",                &SRPDE_Laplacian_2D_Areal::f)
    .method("beta",             &SRPDE_Laplacian_2D_Areal::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_Laplacian_2D_Areal::set_lambda_s)
    .method("set_observations", &SRPDE_Laplacian_2D_Areal::set_observations)
    .method("set_covariates",   &SRPDE_Laplacian_2D_Areal::set_covariates)
    .method("set_subdomains",   &SRPDE_Laplacian_2D_Areal::set_subdomains)
    // solve method
    .method("solve",            &SRPDE_Laplacian_2D_Areal::solve);
}

// Constant Coefficients PDE regularization, sampling at mesh nodes
SRPDE_MODULE(ConstantCoefficients_2D_Order1, Sampling::GeoStatMeshNodes); 
// definition of Rcpp module
typedef R_SRPDE<ConstantCoefficients_2D_Order1, Sampling::GeoStatMeshNodes>
SRPDE_ConstantCoefficients_2D_GeoStatNodes;
RCPP_MODULE(SRPDE_ConstantCoefficients_2D_GeoStatNodes) {
  Rcpp::class_<SRPDE_ConstantCoefficients_2D_GeoStatNodes>("SRPDE_ConstantCoefficients_2D_GeoStatNodes")
    .constructor<ConstantCoefficients_2D_Order1>()
    // getters
    .method("f",                &SRPDE_ConstantCoefficients_2D_GeoStatNodes::f)
    .method("beta",             &SRPDE_ConstantCoefficients_2D_GeoStatNodes::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_ConstantCoefficients_2D_GeoStatNodes::set_lambda_s)
    .method("set_observations", &SRPDE_ConstantCoefficients_2D_GeoStatNodes::set_observations)
    .method("set_covariates",   &SRPDE_ConstantCoefficients_2D_GeoStatNodes::set_covariates)
    // solve method
    .method("solve",            &SRPDE_ConstantCoefficients_2D_GeoStatNodes::solve);
}

// Constant Coefficients PDE regularization, sampling at given locations
SRPDE_MODULE(ConstantCoefficients_2D_Order1, Sampling::GeoStatLocations,
	     void set_locations(const DMatrix<double>& data) { model_.setLocs(data); }
	     )
// definition of Rcpp module
typedef R_SRPDE<ConstantCoefficients_2D_Order1, Sampling::GeoStatLocations>
SRPDE_ConstantCoefficients_2D_GeoStatLocations;
RCPP_MODULE(SRPDE_ConstantCoefficients_2D_GeoStatLocations) {
  Rcpp::class_<SRPDE_ConstantCoefficients_2D_GeoStatLocations>("SRPDE_ConstantCoefficients_2D_GeoStatLocations")
    .constructor<ConstantCoefficients_2D_Order1>()
    // getters
    .method("f",                &SRPDE_ConstantCoefficients_2D_GeoStatLocations::f)
    .method("beta",             &SRPDE_ConstantCoefficients_2D_GeoStatLocations::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_ConstantCoefficients_2D_GeoStatLocations::set_lambda_s)
    .method("set_observations", &SRPDE_ConstantCoefficients_2D_GeoStatLocations::set_observations)
    .method("set_covariates",   &SRPDE_ConstantCoefficients_2D_GeoStatLocations::set_covariates)
    .method("set_locations",    &SRPDE_ConstantCoefficients_2D_GeoStatLocations::set_locations)
    // solve method
    .method("solve",            &SRPDE_ConstantCoefficients_2D_GeoStatLocations::solve);
}

// Constant Coefficients PDE regularization, areal sampling
SRPDE_MODULE(ConstantCoefficients_2D_Order1, Sampling::Areal,
	     void set_subdomains(const DMatrix<int>& data) { model_.setSubdomains(data); }
	     )
// definition of Rcpp module
typedef R_SRPDE<ConstantCoefficients_2D_Order1, Sampling::Areal>
SRPDE_ConstantCoefficients_2D_Areal;
RCPP_MODULE(SRPDE_ConstantCoefficients_2D_Areal) {
  Rcpp::class_<SRPDE_ConstantCoefficients_2D_Areal>("SRPDE_ConstantCoefficients_2D_Areal")
    .constructor<ConstantCoefficients_2D_Order1>()
    // getters
    .method("f",                &SRPDE_ConstantCoefficients_2D_Areal::f)
    .method("beta",             &SRPDE_ConstantCoefficients_2D_Areal::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_ConstantCoefficients_2D_Areal::set_lambda_s)
    .method("set_observations", &SRPDE_ConstantCoefficients_2D_Areal::set_observations)
    .method("set_covariates",   &SRPDE_ConstantCoefficients_2D_Areal::set_covariates)
    .method("set_subdomains",   &SRPDE_ConstantCoefficients_2D_Areal::set_subdomains)
    // solve method
    .method("solve",            &SRPDE_ConstantCoefficients_2D_Areal::solve);
}

// Space Varying PDE regularization, sampling at mesh nodes
SRPDE_MODULE(SpaceVarying_2D_Order1, Sampling::GeoStatMeshNodes); 
// definition of Rcpp module
typedef R_SRPDE<SpaceVarying_2D_Order1, Sampling::GeoStatMeshNodes>
SRPDE_SpaceVarying_2D_GeoStatNodes;
RCPP_MODULE(SRPDE_SpaceVarying_2D_GeoStatNodes) {
  Rcpp::class_<SRPDE_SpaceVarying_2D_GeoStatNodes>("SRPDE_SpaceVarying_2D_GeoStatNodes")
    .constructor<SpaceVarying_2D_Order1>()
    // getters
    .method("f",                &SRPDE_SpaceVarying_2D_GeoStatNodes::f)
    .method("beta",             &SRPDE_SpaceVarying_2D_GeoStatNodes::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_SpaceVarying_2D_GeoStatNodes::set_lambda_s)
    .method("set_observations", &SRPDE_SpaceVarying_2D_GeoStatNodes::set_observations)
    .method("set_covariates",   &SRPDE_SpaceVarying_2D_GeoStatNodes::set_covariates)
    // solve method
    .method("solve",            &SRPDE_SpaceVarying_2D_GeoStatNodes::solve);
}

// Space Varying PDE regularization, sampling at given locations
SRPDE_MODULE(SpaceVarying_2D_Order1, Sampling::GeoStatLocations,
	     void set_locations(const DMatrix<double>& data) { model_.setLocs(data); }
	     )
// definition of Rcpp module
typedef R_SRPDE<SpaceVarying_2D_Order1, Sampling::GeoStatLocations>
SRPDE_SpaceVarying_2D_GeoStatLocations;
RCPP_MODULE(SRPDE_SpaceVarying_2D_GeoStatLocations) {
  Rcpp::class_<SRPDE_SpaceVarying_2D_GeoStatLocations>("SRPDE_SpaceVarying_2D_GeoStatLocations")
    .constructor<SpaceVarying_2D_Order1>()
    // getters
    .method("f",                &SRPDE_SpaceVarying_2D_GeoStatLocations::f)
    .method("beta",             &SRPDE_SpaceVarying_2D_GeoStatLocations::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_SpaceVarying_2D_GeoStatLocations::set_lambda_s)
    .method("set_observations", &SRPDE_SpaceVarying_2D_GeoStatLocations::set_observations)
    .method("set_covariates",   &SRPDE_SpaceVarying_2D_GeoStatLocations::set_covariates)
    .method("set_locations",    &SRPDE_SpaceVarying_2D_GeoStatLocations::set_locations)
    // solve method
    .method("solve",            &SRPDE_SpaceVarying_2D_GeoStatLocations::solve);
}

// Space Varying PDE regularization, areal sampling
SRPDE_MODULE(SpaceVarying_2D_Order1, Sampling::Areal,
	     void set_subdomains(const DMatrix<int>& data) { model_.setSubdomains(data); }
	     )
// definition of Rcpp module
typedef R_SRPDE<SpaceVarying_2D_Order1, Sampling::Areal>
SRPDE_SpaceVarying_2D_Areal;
RCPP_MODULE(SRPDE_SpaceVarying_2D_Areal) {
  Rcpp::class_<SRPDE_SpaceVarying_2D_Areal>("SRPDE_SpaceVarying_2D_Areal")
    .constructor<SpaceVarying_2D_Order1>()
    // getters
    .method("f",                &SRPDE_SpaceVarying_2D_Areal::f)
    .method("beta",             &SRPDE_SpaceVarying_2D_Areal::beta)
    // setters
    .method("set_lambda_s",     &SRPDE_SpaceVarying_2D_Areal::set_lambda_s)
    .method("set_observations", &SRPDE_SpaceVarying_2D_Areal::set_observations)
    .method("set_covariates",   &SRPDE_SpaceVarying_2D_Areal::set_covariates)
    .method("set_subdomains",   &SRPDE_SpaceVarying_2D_Areal::set_subdomains)
    // solve method
    .method("solve",            &SRPDE_SpaceVarying_2D_Areal::solve);
}
