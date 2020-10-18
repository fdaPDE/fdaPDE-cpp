/*
 * param_functors.hpp
 *
 *  Created on: Jun 2, 2015
 *      Author: eardi
 */
#ifndef __PARAM_FUNCTORS_H__
#define __PARAM_FUNCTORS_H__

//#include "matrix_assembler.hpp"

class Function
{
public:
	//virtual ~Function();
	virtual Real operator()(UInt globalNodeIndex, UInt ic = 0) const {return 0;}
};

class Diffusivity
{
	std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > > K_;
public:
	Diffusivity():
			K_(){};
	Diffusivity(const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K):
		K_(K){};
	Diffusivity(SEXP RGlobalVector)
	{
		UInt num_int_nodes = Rf_length(RGlobalVector)/4;
		K_.resize(num_int_nodes);
		for(auto l=0; l<num_int_nodes;l++)
		{
			for(auto j = 0; j < 2; ++j)
			{
				for(auto i = 0; i < 2; ++i)
					K_[l](i,j) = REAL(RGlobalVector)[l*4 + 2*j + i];
			}
		}
	}


	Eigen::Matrix<Real,2,2> operator()(UInt globalNodeIndex, UInt ic1 = 0) const
	{
		return K_[globalNodeIndex];
	}
};

class Advection : public virtual Function
{
	std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > > beta_;
public:
	Advection(const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta):
		beta_(beta){};
	Advection(SEXP RGlobalVector)
	{
		UInt num_int_nodes = Rf_length(RGlobalVector)/2;
		beta_.resize(num_int_nodes);
		for(auto l=0; l<num_int_nodes;l++)
		{
			for(auto j = 0; j < 2; ++j)
			{
					beta_[l](j) = REAL(RGlobalVector)[l*2 + j];
			}
		}
	}

	Real operator() (UInt globalNodeIndex, UInt ic = 0) const
	{
		return beta_[globalNodeIndex][ic];
	}
};

class Reaction : public virtual Function
{
	std::vector<Real> c_;
public:
	Reaction(const std::vector<Real>& c, UInt ic = 0):
		c_(c){};

	Reaction(SEXP RGlobalVector)
	{
		UInt num_int_nodes = Rf_length(RGlobalVector);
		c_.resize(num_int_nodes);
		for(auto l=0; l<num_int_nodes;l++)
		{
					c_[l] = REAL(RGlobalVector)[l];
		}
	}

	Real operator()(UInt globalNodeIndex, UInt ic = 0) const
	{
		return c_[globalNodeIndex];
	}
};

class ForcingTerm : public virtual Function
{
	std::vector<Real> u_;
public:
	ForcingTerm(const std::vector<Real>& u, UInt ic = 0):
		u_(u){};

	ForcingTerm(SEXP RGlobalVector)
	{
		UInt num_int_nodes = Rf_length(RGlobalVector);
		u_.resize(num_int_nodes);
		for(auto l=0; l<num_int_nodes;l++)
		{
					u_[l] = REAL(RGlobalVector)[l];
		}
	}

	Real operator()(UInt globalNodeIndex, UInt ic = 0) const
	{
		return u_[globalNodeIndex];
	}
};


#endif /* PARAM_FUNCTORS_H_ */
