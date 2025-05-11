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

#ifndef __FDAPDE_SOLVERS_MODULE_H__
#define __FDAPDE_SOLVERS_MODULE_H__

// clang-format off

// include core
#include <fdaPDE/core.h>

#include "src/formula.h"
#include "src/distributions.h"
#include "src/solvers/utility.h"

namespace fdapde {

struct ls_solver { };
struct de_solver { };
  
}   // namespace fdapde

// least square solvers
#include "src/solvers/fe_ls_elliptic.h"
#include "src/solvers/fe_ls_separable.h"
#include "src/solvers/fe_ls_parabolic.h"

// density estimation solvers
#include "src/solvers/fe_de_elliptic.h"
#include "src/solvers/fe_de_separable.h"

// clang-format on

#endif   // __FDAPDE_SOLVERS_MODULE_H__
