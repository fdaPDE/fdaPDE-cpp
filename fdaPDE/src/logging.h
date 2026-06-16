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

#ifndef __FDAPDE_LOGGING_H__
#define __FDAPDE_LOGGING_H__

#include <ios>
#include <iostream>
#include <ostream>

namespace fdapde {

#ifdef FDAPDE_ENABLE_COUT

inline std::ostream& cout = std::cout;

#else

class null_ostream {
   public:
    using ostream_manipulator = std::ostream& (*)(std::ostream&);
    using ios_manipulator = std::ios& (*)(std::ios&);
    using ios_base_manipulator = std::ios_base& (*)(std::ios_base&);

    template <typename T>
    constexpr const null_ostream& operator<<(T&&) const noexcept {
        return *this;
    }

    constexpr const null_ostream& operator<<(ostream_manipulator) const noexcept {
        return *this;
    }
    constexpr const null_ostream& operator<<(ios_manipulator) const noexcept {
        return *this;
    }
    constexpr const null_ostream& operator<<(ios_base_manipulator) const noexcept {
        return *this;
    }
};

inline constexpr null_ostream cout {};

#endif

}   // namespace fdapde

#endif   // __FDAPDE_LOGGING_H__
