#ifndef FDAPDE_TEST_IPOPT_OPTIONS_H
#define FDAPDE_TEST_IPOPT_OPTIONS_H

#include <cstdlib>
#include <fstream>

inline void write_ipopt_options() {
    std::ofstream options("ipopt.opt");
    options
        << "print_level 0\n"
        << "sb yes\n"
        << "print_user_options no\n"
        << "print_timing_statistics no\n"
        << "\n";
    if (const char* hsllib = std::getenv("FDAPDE_TEST_HSL_LIBRARY"); hsllib && *hsllib) {
        options
            << "linear_solver ma57\n"
            << "hsllib " << hsllib << "\n\n";
    }
    options
      /*<< "hessian_approximation exact\n"
      << "nlp_scaling_method none\n"
      << "\n"
      << "mu_strategy adaptive\n"
      << "\n"
      << "tol 1e-9\n"
      << "\n"
      << "acceptable_tol 1e-6\n"
      << "acceptable_iter 10\n"
      << "\n"
      << "bound_push 1e-12\n"
      << "bound_frac 1e-12\n"
      << "bound_relax_factor 0\n"*/;
}

#endif // FDAPDE_TEST_IPOPT_OPTIONS_H
