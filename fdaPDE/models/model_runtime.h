// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __MODEL_RUNTIME_H__
#define __MODEL_RUNTIME_H__

namespace fdapde {

// set of possible runtime properties of a model
enum runtime_status {
    // common models flags
    is_lambda_changed,                     // asserted true whenever a new value of \lambda is set
    require_pde_init,                      // if asserted true, pde is not in a valid state and should be reinitialized
    require_penalty_init,                  // if asserted true, the penalty stack must incur a re-initialization
    require_functional_basis_evaluation,   // if asserted true, \Psi matrix must be recomputed
    require_psi_correction,                // if asserted true, matrix \Psi must be updated
    require_data_stack_update,             // asserted true whenever new data are set
    // regression flags
    require_W_update   // asserted true whenever the weights matrix W is modified

};

namespace models {

class model_runtime_handler {
    mutable std::unordered_map<int, bool> runtime_;
   public:
    model_runtime_handler() :
        runtime_ {
	  // we do not force a pde initialization by default since it is already initialized in ModelBase constructor
          {runtime_status::require_penalty_init,                true},
          {runtime_status::require_functional_basis_evaluation, true}
	} {};

    void set(int flag) { runtime_[flag] = true; }
    bool query(int flag) const {
        if (runtime_.find(flag) == runtime_.end()) return false;   // flag was never set
        bool b = runtime_.at(flag);
        runtime_[flag] = false;
        return b;
    }
};

}   // namespace models
}   // namespace fdapde

#endif   // __MODEL_RUNTIME_H__
