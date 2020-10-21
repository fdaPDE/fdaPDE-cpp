/*
 * mesh_object.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */

#include "../Include/Mesh_Objects.h"


const UInt Identifier::NVAL=std::numeric_limits<UInt>::max();


void Point::print(std::ostream & out) const
{
	out<<"Point Id -"<< id_ <<"- "<<"("<<coord_[0]<<","<<coord_[1]<<","<<coord_[2]<<")"<<std::endl<<"------"<<std::endl;
}
void Edge::print(std::ostream & out) const
{
	out<<"Edge Id -"<< id_ <<"- "<<"("<<points_[0].getId()<<","<<points_[1].getId()<<")"<<std::endl;
}
