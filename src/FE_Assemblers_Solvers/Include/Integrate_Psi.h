#ifndef __INTEGRATEPSI_H__
#define __INTEGRATEPSI_H__

#include "../../FdaPDE.h"

template <UInt Nodes, UInt mydim, UInt ndim>
inline Real integratePsi(const Element<Nodes,mydim,ndim>& t, UInt k)
{
	std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
	return 0;
}

template<>
inline Real integratePsi(const Element<3,2,2>& t, UInt k)
{
	return t.getArea()/3;
}

template<>
inline Real integratePsi(const Element<6,2,2>& t, UInt k)
{
	if (k==3 || k==4 || k==5)
		return t.getArea()/3;
	else
		return 0;
}

template<>
inline Real integratePsi(const Element<3,2,3>& t, UInt k)
{
	return t.getArea()/3;
}

template<>
inline Real integratePsi(const Element<6,2,3>& t, UInt k)
{
	if (k==3 || k==4 || k==5)
		return t.getArea()/3;
	else
		return 0;
}

template<>
inline Real integratePsi(const Element<4,3,3>& t, UInt k)
{
	return t.getVolume()/4;
}

#endif
