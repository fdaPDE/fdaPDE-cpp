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

#ifndef __FORMULA_H__
#define __FORMULA_H__

namespace fdapde {

class Formula {
    std::string lhs_;
    std::vector<std::string> rhs_;

    void throw_parse_error_(const std::string& msg) {
        throw std::runtime_error(std::string("Formula parse error: ") + msg);
    }
   public:
    Formula() noexcept = default;
    Formula(const std::string& formula) {
        std::size_t tilde = formula.find('~');
        if (tilde == std::string::npos) { throw_parse_error_("no '~' delimiter found."); }
        lhs_ = formula.substr(0, tilde);
        std::erase(lhs_, ' ');
        if (lhs_.empty()) { throw_parse_error_("no lhs found."); }
        // rhs parsing logic
        std::string rhs = formula.substr(tilde + 1, formula.size() - tilde - 1);
        std::erase(rhs, ' ');
	if (rhs.empty())  { throw_parse_error_("no rhs found."); }
        size_t pos = rhs.find('+');
        std::string token;
        while (pos != std::string::npos) {
            token = rhs.substr(0, pos);
            rhs_.push_back(token);
            rhs.erase(0, pos + 1);
	    pos = rhs.find('+');
        }
        rhs_.push_back(rhs);
    }
    // observers
    const std::string& lhs() const { return lhs_; }
    const std::vector<std::string>& rhs() const { return rhs_; }
};

}   // namespace fdapde

#endif // __FORMULA_H__
