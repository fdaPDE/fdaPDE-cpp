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
    // available tokens
    struct efx_token {
       private:
        std::string cov_;
        std::string efx_;
       public:
        efx_token(const std::string& cov, const std::string& efx) : cov_(cov), efx_(efx) { }
        // observers
        const std::string& cov() const { return cov_; }
        const std::string& efx() const { return efx_; }
        friend std::ostream& operator<<(std::ostream& os, const efx_token& m) {
            os << "(" << m.cov_ << ", " << m.efx_ << ")";
            return os;
        }
    };

    std::string lhs_;
    std::vector<std::string> covs_;
    std::vector<efx_token> efxs_;

    void throw_parse_error_(const std::string& msg) {
        throw std::runtime_error(std::string("Formula parse error: ") + msg);
    }
    void analyze_token_(std::string token) {
        size_t pos = token.find('|');
        if (pos != std::string::npos) {
            std::string cov_token, efx_token;
            cov_token = token.substr(0, pos);
            token.erase(0, pos + 1);
            efx_token = token;
            efxs_.emplace_back(cov_token, efx_token);
        } else {
            covs_.push_back(token);
        }
        return;
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
        if (rhs.empty()) { throw_parse_error_("no rhs found."); }
        size_t pos = rhs.find('+');
        std::string token;
        while (pos != std::string::npos) {
            token = rhs.substr(0, pos);
            analyze_token_(token);
            rhs.erase(0, pos + 1);
            pos = rhs.find('+');
        }
        analyze_token_(rhs);
    }
    // observers
    const std::string& lhs() const { return lhs_; }
    const std::vector<std::string>& covs() const { return covs_; }
    const std::vector<efx_token>& efxs() const { return efxs_; }
};

}   // namespace fdapde

#endif // __FORMULA_H__
