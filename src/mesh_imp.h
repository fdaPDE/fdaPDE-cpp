#ifndef MESH_IMP_H_
#define MESH_IMP_H_

#include<iostream>
#include<fstream>
#include<sstream>

#ifdef R_VERSION_
template <UInt ORDER>
MeshHandler<ORDER,2,2>::MeshHandler(SEXP mesh, UInt search)
{
	mesh_ 		= mesh;
	points_ 	= REAL(VECTOR_ELT(mesh_, 0));
	elements_  = INTEGER(VECTOR_ELT(mesh_, 3));
	edges_ 		= INTEGER(VECTOR_ELT(mesh_, 6));
	neighbors_  = INTEGER(VECTOR_ELT(mesh_, 8));

	num_nodes_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 0), R_DimSymbol))[0];
	num_elements_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 3), R_DimSymbol))[0];
	num_edges_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 6), R_DimSymbol))[0];
	search_ = search;
	
	if (search == 2) { //if tree search, construct a tree mesh
		// Rprintf("mesh TYPE: %d \n",TYPEOF(mesh_)); //VECSXP, list (generic vector), 19
		// Rprintf("mesh LENGTH: %d \n",XLENGTH(mesh_));
		
		int mesh_len = XLENGTH(mesh_);
		if (mesh_len == 11) { //don't have tree mesh information (length==11)
			ADTree<Element<3*ORDER,2,2>> tmp(points_, elements_, num_nodes_, num_elements_);
			tree_ = tmp;
		} else {
			//RECIEVE TREE INFORMATION FROM R
			//tree_header information
			int tree_loc_ = num_elements_;
			int tree_lev_ = INTEGER(VECTOR_ELT(mesh_, 11))[0];
			int ndimp_ = 2;
			int ndimt_ = 4;
			int nele_ = num_elements_;
			int iava_ = num_elements_+1;
			int iend_ = num_elements_+1;

			std::vector<Real>  origin_;
			origin_.assign(REAL(VECTOR_ELT(mesh_, 12)), REAL(VECTOR_ELT(mesh_, 12))+ndimt_);
			std::vector<Real> scalingfactors_;
			scalingfactors_.assign(REAL(VECTOR_ELT(mesh_, 13)), REAL(VECTOR_ELT(mesh_, 13))+ndimt_);
			
			Domain<Element<3*ORDER, 2, 2>> tree_domain(origin_, scalingfactors_);
			TreeHeader<Element<3*ORDER,2, 2>> tree_header(tree_loc_, tree_lev_, ndimp_, ndimt_, nele_, iava_, iend_, tree_domain);
			

			//treenode information (number of nodes = number of elements+1)
			std::vector<Id> id_;
			id_.assign(INTEGER(VECTOR_ELT(mesh_, 14)), INTEGER(VECTOR_ELT(mesh_, 14))+num_elements_+1);
			std::vector<int> node_left_child_;
			node_left_child_.assign(INTEGER(VECTOR_ELT(mesh_, 15)), INTEGER(VECTOR_ELT(mesh_, 15))+num_elements_+1);
			std::vector<int> node_right_child_;
			node_right_child_.assign(INTEGER(VECTOR_ELT(mesh_, 16)), INTEGER(VECTOR_ELT(mesh_, 16))+num_elements_+1);
			Real* box_ = REAL(VECTOR_ELT(mesh_, 17));

			UInt num_tree_nodes = id_.size();
			std::vector<TreeNode<Element<3*ORDER,2,2>>> tree_nodes;
			for (UInt i=0; i<num_tree_nodes; i++) {
				std::vector<Real> coord;
				for (UInt j=0; j<ndimt_; j++) {
					coord.push_back(box_[i + num_tree_nodes*j]);
				}
				Box<2> box (coord);
				TreeNode<Element<3*ORDER,2,2>> tree_node(box, id_[i], node_left_child_[i], node_right_child_[i]);
				tree_nodes.push_back(tree_node);
			}
			

			ADTree<Element<3*ORDER,2,2>> tmp(tree_header, tree_nodes);
			tree_ = tmp;
		}
	} //if tree search, construct a tree mesh
}
#endif

template <UInt ORDER>
Point MeshHandler<ORDER,2,2>::getPoint(Id id)
{
	Point point(id, Identifier::NVAL, points_[id], points_[num_nodes_+id]);
	return point;
}

template <UInt ORDER>
Edge MeshHandler<ORDER,2,2>::getEdge(Id id)
{
	Id id_start_point = edges_[id];
	Id id_end_point = edges_[num_edges_+id];
	Edge edge(id, Identifier::NVAL, Point(id_start_point, Identifier::NVAL, points_[id_start_point], points_[num_nodes_+id_start_point]),
						Point(id_end_point, Identifier::NVAL, points_[id_end_point], points_[num_nodes_+id_end_point]));
	return edge;
}

template <UInt ORDER>
Element<3*ORDER,2,2> MeshHandler<ORDER,2,2>::getElement(Id id) const
{
	std::vector<Point> element_points;
	element_points.resize(3*ORDER);
	Id id_current_point;
	for (int i=0; i<3*ORDER; ++i)
	{
		id_current_point = elements_[i*num_elements_ + id];
		element_points[i]= Point(id_current_point, 
								 Identifier::NVAL, 
								 points_[id_current_point],
								 points_[num_nodes_+id_current_point]);
						  
								  
												
	}
	return Element<3*ORDER,2,2>(id, element_points);
}

template <UInt ORDER>
Element<3*ORDER,2,2>MeshHandler<ORDER,2,2>::getNeighbors(Id id_element, UInt number) const
{
	Id id_neighbour = neighbors_[number * num_elements_ + id_element];
	//std::cout<<"Neighbour id "<< id_neighbour;
	if (id_neighbour == -1) return Element<3*ORDER,2,2>(); //Element with NVAL ID

	return getElement(id_neighbour);
}

template <UInt ORDER>
Element<3*ORDER,2,2> MeshHandler<ORDER,2,2>::findLocationNaive(Point point) const
{
	Element<3*ORDER,2,2> current_element;
	//std::cout<<"Start searching point naively \n";
	for(Id id=0; id < num_elements_; ++id)
	{
		current_element = getElement(id);
		if(current_element.isPointInside(point))
			return current_element;
	}
	//std::cout<<"Point not found \n";
	return Element<3*ORDER,2,2>(); //default element with NVAL ID
}

// Visibility walk algorithm which uses barycentric coordinate [Sundareswara et al]
//Starting triangles usually n^(1/3) points
template <UInt ORDER>
Element<3*ORDER,2,2> MeshHandler<ORDER,2,2>::findLocationWalking(const Point& point, const Element<3*ORDER,2,2>& starting_element) const
{

	//Walking algorithm to the point
	Element<3*ORDER,2,2> current_element = starting_element;

	int direction=0;

	//Test for found Element, or out of border
	while(current_element.getId() != Identifier::NVAL && !current_element.isPointInside(point) )
	{
		direction = current_element.getPointDirection(point);
		//std::cout<<"Direction "<<direction<<";";
		current_element = getNeighbors(current_element.getId(), direction);
  	    //std::cout<<" ID "<<current_element.getId();
	}

	return current_element;
}


template <UInt ORDER>
Element<3*ORDER,2,2> MeshHandler<ORDER,2,2>::findLocationTree(const Point& point) const {
		std::vector<Real> region(4);
		bool result;
		std::set<int> found;
		int index;
		Element<3*ORDER,2,2> tmp;
		region[0] = point[0];	
		region[1] = point[1];
		region[2] = point[0];
		region[3] = point[1];
	
		result = tree_.search(region, found);
		if(result == 0) {
			return Element<3*ORDER,2,2>();
		}

		for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
			index = *i;
			index = this -> tree_.pointId(index);
	  		tmp = this -> getElement(index);
			result = tmp.isPointInside(point);
			if(result == 1) {
				return tmp;
			}
		}
	return Element<3*ORDER,2,2>();
}


template <UInt ORDER>
Real MeshHandler<ORDER,2,2>::elementMeasure(Id id) const
{
	std::vector<Point> p;
	p.resize(3*ORDER);
	Id id_current_point;
	for (int i=0; i<3*ORDER; ++i)
	{
		id_current_point = elements_[i*num_elements_ + id];
		p[i]= Point(id_current_point, Identifier::NVAL, points_[id_current_point],points_[num_nodes_+id_current_point]);
	}
	Real area = std::abs((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[2][0]-p[0][0])*(p[1][1]-p[0][1]))/2;
	return area;
}


/*std::ostream & operator<<(std::ostream & out, MeshHandler const& m){
	out<< " ***** MESH  INFORMATION ******"<<std::endl;
	out<<" Num Points="<<m.num_nodes()<<" "<<" Num elements="<<m.num_triangles()<<" "
			<<"Num. edges="<<m.num_edges();
			//<<" "<<"Num Boundary Edges="<<m.num_bEdges()<<std::endl;
	out<< "POINTS:"<<std::endl;
	int oprec=out.precision(10);
	std::ios_base::fmtflags oflags=
			out.setf(std::ios_base::scientific,std::ios_base::floatfield);
	for (UInt i=0;i<m.num_nodes();++i){
		Point p=m.point(i);
		double x=p[0];
		double y=p[1];
		out<<i<<" "<<x<<" "<<y<<std::endl;
	}
	out<<" TRIANGLE CONNECTIVITY AND AREA:"<<std::endl;
	for (UInt i=0; i<m.num_elements();++i){
		Triangle t=m.triangle(i);
		out<<i<<" "<<t[0].id()<<" "<<t[1].id()<<" "<<t[2].id()<<
		  " "<<t.measure()<<"  Edge: "<<t.getEdges_id(0)<<" "<<t.getEdges_id(1)<<" "<<t.getEdges_id(2)<<std::endl;
	}
	out.precision(oprec);
	out.flags(oflags);
	return out;
}*/

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printPoints(std::ostream & out)
{
	for(UInt i = 0; i < num_nodes_; ++i)
	{
		out<<"-"<< i <<"-"<<"("<<points_[i]<<","<<points_[num_nodes_+i]<<")"<<std::endl<<"------"<<std::endl;
	}
}

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printEdges(std::ostream & out)
{

	out << "Numero lati: "<< num_edges_ <<std::endl;
	for (UInt i = 0; i < num_edges_; ++i )
	{
		out<<"Lato ("<<edges_[i]<<","<<edges_[num_edges_+i]<<")"<<std::endl;
	}

}

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printElements(std::ostream & out)
{

	out << "# Triangles: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < 3*ORDER; ++k)
			out<<elements_[k*num_elements_ + i]<<"   ";
		out<<std::endl;
	}

}

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printNeighbors(std::ostream & out)
{

	out << "# Neighbors list: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < 3*ORDER; ++k)
			out<<neighbors_[k*num_elements_ + i]<<"   ";
		out<<std::endl;
	}

}


template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printTree(std::ostream & out)
{
	
	out << "# Tree characteristic: " <<std::endl;
	out << tree_ << std::endl;
	
}

//////////////////////////////////////////////////////////
// Implementation of class MeshHandler for surface mesh //
//////////////////////////////////////////////////////////


#ifdef R_VERSION_
template <UInt ORDER>
MeshHandler<ORDER,2,3>::MeshHandler(SEXP mesh, UInt search)
{

	mesh_ = mesh;
	num_nodes_ = INTEGER(VECTOR_ELT(mesh_,0))[0];
	num_elements_ = INTEGER(VECTOR_ELT(mesh_,1))[0];
	points_ = REAL(VECTOR_ELT(mesh_, 2));
	elements_ = INTEGER(VECTOR_ELT(mesh_, 3));
	search_ = search;			  
	
	if (search == 2) { //if tree search, construct a tree mesh	
		//Rprintf("mesh LENGTH: %d \n", XLENGTH(mesh_));
		int mesh_len = XLENGTH(mesh_);
		if (mesh_len == 5) { //don't have tree mesh information (length==5)
			ADTree<Element<3*ORDER,2,3>> tmp(points_, elements_, num_nodes_, num_elements_);
			tree_ = tmp;
		} else {
			//RECIEVE TREE INFORMATION FROM R
			//tree_header information
			int tree_loc_ = num_elements_;
			int tree_lev_ = INTEGER(VECTOR_ELT(mesh_, 5))[0];
			int ndimp_ = 3;
			int ndimt_ = 6;
			int nele_ = num_elements_;
			int iava_ = num_elements_+1;
			int iend_ = num_elements_+1;

			std::vector<Real>  origin_;
			origin_.assign(REAL(VECTOR_ELT(mesh_, 6)), REAL(VECTOR_ELT(mesh_, 6))+ndimt_);
			std::vector<Real> scalingfactors_;
			scalingfactors_.assign(REAL(VECTOR_ELT(mesh_, 7)), REAL(VECTOR_ELT(mesh_, 7))+ndimt_);
			
			Domain<Element<3*ORDER, 2, 3>> tree_domain(origin_, scalingfactors_);
			TreeHeader<Element<3*ORDER,2, 3>> tree_header(tree_loc_, tree_lev_, ndimp_, ndimt_, nele_, iava_, iend_, tree_domain);
			

			//treenode information (number of nodes = number of elements+1)
			std::vector<Id> id_;
			id_.assign(INTEGER(VECTOR_ELT(mesh_, 8)), INTEGER(VECTOR_ELT(mesh_, 8))+num_elements_+1);
			std::vector<int> node_left_child_;
			node_left_child_.assign(INTEGER(VECTOR_ELT(mesh_, 9)), INTEGER(VECTOR_ELT(mesh_, 9))+num_elements_+1);
			std::vector<int> node_right_child_;
			node_right_child_.assign(INTEGER(VECTOR_ELT(mesh_, 10)), INTEGER(VECTOR_ELT(mesh_, 10))+num_elements_+1);
			Real* box_ = REAL(VECTOR_ELT(mesh_, 11));

			UInt num_tree_nodes = id_.size();
			std::vector<TreeNode<Element<3*ORDER,2,3>>> tree_nodes;
			for (UInt i=0; i<num_tree_nodes; i++) {
				std::vector<Real> coord;
				for (UInt j=0; j<ndimt_; j++) {
					coord.push_back(box_[i + num_tree_nodes*j]);
				}
				Box<3> box (coord);
				TreeNode<Element<3*ORDER,2,3>> tree_node(box, id_[i], node_left_child_[i], node_right_child_[i]);
				tree_nodes.push_back(tree_node);
			}
			
			ADTree<Element<3*ORDER,2,3>> tmp(tree_header, tree_nodes);
			tree_ = tmp;
		}
	}

}
#endif


template <UInt ORDER>
void MeshHandler<ORDER,2,3>::importfromCSV(std::string &filename){

	UInt nnodes;
	UInt ntriangles;
	UInt point_index;
	std::string line;
	std::string dummy;
	char comma;

	std::ifstream file;
	file.open(filename);


	// Read the number of points
	getline(file,line);
	std::istringstream ss(line);

	ss >> dummy; // throw away "num_points"
	ss >> nnodes;

	num_nodes_ = nnodes;
	// points_.resize(3*nnodes);

	// Read the number of points
	getline(file,line);
	std::istringstream ss2(line);

	ss2 >> dummy; // throw away "num_triangles"
	ss2 >> ntriangles;

	num_elements_ = ntriangles;
	// elements_.resize(3*ORDER*ntriangles);


	getline(file,line); //skip a white line

	// READ THE VERTICES MATRIX
	for(UInt i=0; i<nnodes; ++i){
		std::getline(file,line);
		std::istringstream ss(line);
		ss>>points_[3*i];
		ss>>comma;
		ss>>points_[3*i+1];
		ss>>comma;
		ss>>points_[3*i+2];
	};

	getline(file,line); //skip a white line

	// READ THE CONNECTIVIY MATRIX


	for(UInt i=0; i<ntriangles; ++i){
		std::getline(file,line);
		std::istringstream ss(line);
		ss>>point_index;
		elements_[i*3] = --point_index;
		ss>>comma;
		ss>>point_index;
		elements_[i*3+1] = --point_index;
		ss>>comma;
		ss>>point_index;
		elements_[i*3+2] = --point_index;

	};


};


template <UInt ORDER>
Point MeshHandler<ORDER,2,3>::getPoint(Id id)
{
	Point point(id, Identifier::NVAL, points_[3*id], points_[3*id+1], points_[3*id+2]);
	return point;
}

template <UInt ORDER>
Element<3*ORDER,2,3> MeshHandler<ORDER,2,3>::getElement(Id id) const
{
	std::vector<Point> element_points;
	element_points.resize(3*ORDER);
	Id id_current_point;
	for (int i=0; i<3*ORDER; ++i)
	{
		id_current_point = elements_[3*ORDER * id + i];
		element_points[i]= Point(id_current_point, 
								 Identifier::NVAL, 
								 points_[3*id_current_point],
								 points_[3*id_current_point+1],
								 points_[3*id_current_point+2]);									   
	}
	return Element<3*ORDER,2,3>(id, element_points);	 													
}


template <UInt ORDER>
Element<3*ORDER,2,3> MeshHandler<ORDER,2,3>::findLocationNaive(Point point) const
{
	Element<3*ORDER,2,3> current_element;
	//std::cout<<"Start searching point naively \n";
	for(Id id=0; id < num_elements_; ++id)
	{
		current_element = getElement(id);
		if(current_element.isPointInside(point))
			return current_element;
	}
	//std::cout<<"Point not found \n";
	return Element<3*ORDER,2,3>(); //default element with NVAL ID
}

template <UInt ORDER>
Element<3*ORDER,2,3> MeshHandler<ORDER,2,3>::findLocationTree(const Point& point) const {
	std::vector<Real> region(6);
	bool result;
	std::set<int> found;
	int index;
	Element<3*ORDER,2,3> tmp;
	region[0] = point[0];	
	region[1] = point[1];
	region[2] = point[2];
	region[3] = point[0];
	region[4] = point[1];
	region[5] =	point[2];

	result = tree_.search(region, found);
	if(result == 0) {
		return Element<3*ORDER,2,3>();
	}
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
		index = *i;
		index = this -> tree_.pointId(index);
  		tmp = this -> getElement(index);
		result = tmp.isPointInside(point);
		if(result == 1) {
			return tmp;
		}
	}
	return Element<3*ORDER,2,3>();
}

template <UInt ORDER>
Real MeshHandler<ORDER,2,3>::elementMeasure(Id id) const
{
	std::vector<Point> p;
	p.resize(3*ORDER);
	Id id_current_point;
	for (int i=0; i<3*ORDER; ++i)
	{
		id_current_point = elements_[3*ORDER * id + i];
		p[i]= Point(id_current_point, Identifier::NVAL, points_[3*id_current_point],points_[3*id_current_point+1],points_[3*id_current_point+2]);
					   
								 
								   
									
	}
	Real a2 = std::pow(p[1][0]-p[2][0],2)+std::pow(p[1][1]-p[2][1],2)+std::pow(p[1][2]-p[2][2],2);
	Real b2 = std::pow(p[0][0]-p[2][0],2)+std::pow(p[0][1]-p[2][1],2)+std::pow(p[0][2]-p[2][2],2);
	Real c2 = std::pow(p[0][0]-p[1][0],2)+std::pow(p[0][1]-p[1][1],2)+std::pow(p[0][2]-p[1][2],2);
	Real area = std::sqrt(4*(a2*b2+a2*c2+b2*c2)-std::pow(a2+b2+c2,2))/4; //Heron's formula
	return area;
}

template <UInt ORDER>
void MeshHandler<ORDER,2,3>::printPoints(std::ostream & out)
{
std::cout<<"printing points"<<"\n";
	for(UInt i = 0; i < num_nodes_; ++i)
	{
		out<<"-"<< i <<"-"<<"("<<points_[3*i]<<","<<points_[3*i+1]<<","<<points_[3*i+2]<<")"<<std::endl<<"------"<<std::endl;
	}
}

template <UInt ORDER>
void MeshHandler<ORDER,2,3>::printElements(std::ostream & out)
{

	out << "# Triangles: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < ORDER * 3; ++k)
			out<<elements_[i*3*ORDER + k]<<"   ";
		out<<std::endl;
	}

}



//////////////////////////////////////////////////////////
// Implementation of class MeshHandler for volume mesh //
//////////////////////////////////////////////////////////

#ifdef R_VERSION_
template <UInt ORDER>
MeshHandler<ORDER,3,3>::MeshHandler(SEXP mesh, UInt search)
{

	mesh_ = mesh;
	num_nodes_ = INTEGER(VECTOR_ELT(mesh_,0))[0];
	num_elements_ = INTEGER(VECTOR_ELT(mesh_,1))[0];
	points_ = REAL(VECTOR_ELT(mesh_, 2));
	elements_ = INTEGER(VECTOR_ELT(mesh_, 3));
	search_ = search;
	
	if (search == 2) { //if tree search, construct a tree mesh
		// Rprintf("mesh LENGTH: %d \n",XLENGTH(mesh_));
		int mesh_len = XLENGTH(mesh_);
		if (mesh_len == 5) { //don't have tree mesh information (length==5)
			ADTree<Element<6*ORDER-2,3,3>> tmp(points_, elements_, num_nodes_, num_elements_);
			tree_ = tmp;
		} else {
			//RECIEVE TREE INFORMATION FROM R
			//tree_header information
			int tree_loc_ = num_elements_;
			int tree_lev_ = INTEGER(VECTOR_ELT(mesh_, 5))[0];
			int ndimp_ = 3;
			int ndimt_ = 6;
			int nele_ = num_elements_;
			int iava_ = num_elements_+1;
			int iend_ = num_elements_+1;

			std::vector<Real>  origin_;
			origin_.assign(REAL(VECTOR_ELT(mesh_, 6)), REAL(VECTOR_ELT(mesh_, 6))+ndimt_);
			std::vector<Real> scalingfactors_;
			scalingfactors_.assign(REAL(VECTOR_ELT(mesh_, 7)), REAL(VECTOR_ELT(mesh_, 7))+ndimt_);
			
			Domain<Element<6*ORDER-2,3,3>> tree_domain(origin_, scalingfactors_);
			TreeHeader<Element<6*ORDER-2,3,3>> tree_header(tree_loc_, tree_lev_, ndimp_, ndimt_, nele_, iava_, iend_, tree_domain);
			

			//treenode information (number of nodes = number of elements+1)
			std::vector<Id> id_;
			id_.assign(INTEGER(VECTOR_ELT(mesh_, 8)), INTEGER(VECTOR_ELT(mesh_, 8))+num_elements_+1);
			std::vector<int> node_left_child_;
			node_left_child_.assign(INTEGER(VECTOR_ELT(mesh_, 9)), INTEGER(VECTOR_ELT(mesh_, 9))+num_elements_+1);
			std::vector<int> node_right_child_;
			node_right_child_.assign(INTEGER(VECTOR_ELT(mesh_, 10)), INTEGER(VECTOR_ELT(mesh_, 10))+num_elements_+1);
			Real* box_ = REAL(VECTOR_ELT(mesh_, 11));

			UInt num_tree_nodes = id_.size();
			std::vector<TreeNode<Element<6*ORDER-2,3,3>>> tree_nodes;
			for (UInt i=0; i<num_tree_nodes; i++) {
				std::vector<Real> coord;
				for (UInt j=0; j<ndimt_; j++) {
					coord.push_back(box_[i + num_tree_nodes*j]);
				}
				Box<3> box (coord);
				TreeNode<Element<6*ORDER-2,3,3>> tree_node(box, id_[i], node_left_child_[i], node_right_child_[i]);
				tree_nodes.push_back(tree_node);
			}
			
			ADTree<Element<6*ORDER-2,3,3>> tmp(tree_header, tree_nodes);
			tree_ = tmp;
		}
	}
}
#endif


template <UInt ORDER>
Point MeshHandler<ORDER,3,3>::getPoint(Id id)
{
	Point point(id, Identifier::NVAL, points_[3*id], points_[3*id+1], points_[3*id+2]);
	return point;
}

template <UInt ORDER>
Element<6*ORDER-2,3,3> MeshHandler<ORDER,3,3>::getElement(Id id) const
{
	std::vector<Point> element_points;
	element_points.resize(6*ORDER-2);
	Id id_current_point;
	for (int i=0; i<6*ORDER-2; ++i)
	{
		id_current_point = elements_[(6*ORDER-2) * id + i];
		element_points[i]= Point(id_current_point, 
								 Identifier::NVAL, 
								 points_[3*id_current_point],
								 points_[3*id_current_point+1],
								 points_[3*id_current_point+2]);									   
	}
	return Element<6*ORDER-2,3,3>(id, element_points);													  
}

template <UInt ORDER>
Element<6*ORDER-2,3,3> MeshHandler<ORDER,3,3>::findLocationNaive(Point point) const
{
	Element<6*ORDER-2,3,3> current_element;
	//std::cout<<"Start searching point naively \n";
	for(Id id=0; id < num_elements_; ++id)
	{
		current_element = getElement(id);
		if(current_element.isPointInside(point))
			return current_element;
	}
	//std::cout<<"Point not found \n";
	return Element<6*ORDER-2,3,3>(); //default element with NVAL ID
}

template <UInt ORDER>
Element<6*ORDER-2,3,3> MeshHandler<ORDER,3,3>::findLocationTree(const Point& point) const {
	std::vector<Real> region(6);
	bool result;
	std::set<int> found;
	int index;
	Element<6*ORDER-2,3,3> tmp;
	region[0] = point[0];	
	region[1] = point[1];
	region[2] = point[2];
	region[3] = point[0];
	region[4] = point[1];
	region[5] =	point[2];

	result = tree_.search(region, found);
	if(result == 0) {
		return Element<6*ORDER-2,3,3>();
	}
	for (std::set<int>::iterator i = found.begin(); i != found.end(); i++) {
		index = *i;
		index = this -> tree_.pointId(index);
  		tmp = this -> getElement(index);
		result = tmp.isPointInside(point);
		if(result == 1) {
			return tmp;
		}
	}
	return Element<6*ORDER-2,3,3>();
}

template <UInt ORDER>
Real MeshHandler<ORDER,3,3>::elementMeasure(Id id) const
{
	std::vector<Point> p;
	p.resize(6*ORDER-2);
	Id id_current_point;
	for (int i=0; i<6*ORDER-2; ++i)
	{
		id_current_point = elements_[(6*ORDER-2) * id + i];
		p[i]= Point(id_current_point, Identifier::NVAL, points_[3*id_current_point],points_[3*id_current_point+1],points_[3*id_current_point+2]);
					   
								 
								   
									
	}
	Real volume = std::abs((p[1][0]-p[0][0])*((p[2][1]-p[0][1])*(p[3][2]-p[0][2])-(p[3][1]-p[0][1])*(p[2][2]-p[0][2]))-(p[2][0]-p[0][0])*((p[1][1]-p[0][1])*(p[3][2]-p[0][2])-(p[3][1]-p[0][1])*(p[1][2]-p[0][2]))+(p[3][0]-p[0][0])*((p[1][1]-p[0][1])*(p[2][2]-p[0][2])-(p[2][1]-p[0][1])*(p[1][2]-p[0][2])))/6;
	return volume;
}

template <UInt ORDER>
void MeshHandler<ORDER,3,3>::printPoints(std::ostream & out)
{
std::cout<<"printing points"<<"\n";
	for(UInt i = 0; i < num_nodes_; ++i)
	{
		out<<"-"<< i <<"-"<<"("<<points_[3*i]<<","<<points_[3*i+1]<<","<<points_[3*i+2]<<")"<<std::endl<<"------"<<std::endl;
	}
}

template <UInt ORDER>
void MeshHandler<ORDER,3,3>::printElements(std::ostream & out)
{

	out << "# Tetrahedrons: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < (6*ORDER-2); ++k)
			out<<elements_[i*(6*ORDER-2) + k]<<"   ";
		out<<std::endl;
	}

}

#endif