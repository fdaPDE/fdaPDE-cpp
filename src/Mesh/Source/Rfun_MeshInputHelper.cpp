#include "../../FdaPDE.h"

#include "../Include/Mesh_Input_Helper.h"

#include <array>

extern "C" {


SEXP CPP_SurfaceMeshHelper(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};


  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 4));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    edges_list.assemble_subs(result, 0);
    edges_list.mark_boundary(result, 1);
    edges_list.compute_neighbors(result, 3);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);

	UNPROTECT(1);

  return result;
}


SEXP CPP_SurfaceMeshOrder2(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 6));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    edges_list.assemble_subs(result, 0);
    edges_list.mark_boundary(result, 1);
    edges_list.compute_neighbors(result, 3);
    edges_list.order2extend(result, 5);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);
  compute_midpoints(result, Rnodes, 4, 0);

	UNPROTECT(1);

  return result;
}

SEXP CPP_TriangleMeshSplit(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 2));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    split(result, Rtriangles, 0, edges_list);
    compute_midpoints(result, Rnodes, 1, edges_list);
  }

	UNPROTECT(1);

  return result;
}

SEXP CPP_TriangleMeshSplitOrder2(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 1));

  {
    simplex_container<2> edges_list(Rtriangles, Rnodes, EDGES_ORDERING);
    split(result, Rtriangles, 0, edges_list);
  }

  UNPROTECT(1);

  return result;
}



SEXP CPP_VolumeMeshHelper(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> FACES_ORDERING = {1,2,3,0,2,3,0,1,3,0,1,2};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 4));

  {
    simplex_container<3> faces_list(Rtetrahedrons, Rnodes, FACES_ORDERING);
    faces_list.assemble_subs(result, 0);
    faces_list.mark_boundary(result, 1);
    faces_list.compute_neighbors(result, 3);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);

	UNPROTECT(1);

  return result;

}

SEXP CPP_VolumeMeshOrder2(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> FACES_ORDERING = {1,2,3,0,2,3,0,1,3,0,1,2};
  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 6));

  {
    simplex_container<3> faces_list(Rtetrahedrons, Rnodes, FACES_ORDERING);
    faces_list.assemble_subs(result, 0);
    faces_list.mark_boundary(result, 1);
    faces_list.compute_neighbors(result, 3);
  }

  mark_boundary_nodes(result, Rnodes, 2, 0, 1);

  {
    simplex_container<2> edges_list(Rtetrahedrons, Rnodes, EDGES_ORDERING);
    edges_list.order2extend(result, 5);
    compute_midpoints(result, Rnodes, 4, edges_list);
  }

	UNPROTECT(1);

  return result;

}


SEXP CPP_TetraMeshSplit(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 2));

  {
    simplex_container<2> edges_list(Rtetrahedrons, Rnodes, EDGES_ORDERING);
    split3D(result, Rtetrahedrons, 0, edges_list);
    compute_midpoints(result, Rnodes, 1, edges_list);
  }

	UNPROTECT(1);

  return result;
}

SEXP CPP_TetraMeshSplitOrder2(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};

  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 1));

  {
    simplex_container<2> edges_list(Rtetrahedrons, Rnodes, EDGES_ORDERING);
    split3D(result, Rtetrahedrons, 0, edges_list);
  }

  UNPROTECT(1);

  return result;
}

SEXP CPP_EdgeMeshHelper(SEXP Redges, SEXP Rnodes){
    
    static constexpr std::array<UInt, 2> NODES_ORDERING = {1,0};
    
    
    SEXP result = NILSXP;
    result = PROTECT(Rf_allocVector(VECSXP, 4));
    
    {
      simplex_container<1> nodes_list(Redges, Rnodes, NODES_ORDERING);
      nodes_list.assemble_subs(result, 0);
      nodes_list.mark_boundary(result, 1);
      // compute_neighbors fills the index-th and the (index+1)-th positions in resultf
      nodes_list.compute_neighbors(result, 2);
    }
    
    UNPROTECT(1);
    
    return result;
}

SEXP CPP_EdgeMeshOrder2(SEXP Redges, SEXP Rnodes){

    static constexpr std::array<UInt, 2> NODES_ORDERING = {1,0};


    SEXP result = NILSXP;
    result = PROTECT(Rf_allocVector(VECSXP, 6));

    {
        simplex_container<1> nodes_list(Redges, Rnodes, NODES_ORDERING);
        nodes_list.assemble_subs(result, 0);
        nodes_list.mark_boundary(result, 1);
        //compute_neighbors fills index 2 and index 3
        nodes_list.compute_neighbors(result, 2);
        compute_midpoints(result, Rnodes,Redges,4);

        //midpoints global numbering
        SET_VECTOR_ELT(result,5, Rf_allocMatrix(INTSXP, nodes_list.get_num_elements(),1 ));
        RIntegerMatrix midpoints(VECTOR_ELT(result,5));
        UInt num_points = nodes_list.get_num_points();
        for(UInt i=0; i<nodes_list.get_num_elements(); ++i, ++num_points)
            midpoints[i]=num_points;
    }


    UNPROTECT(1);

    return result;
}

SEXP CPP_EdgeMeshSplit(SEXP Redges, SEXP Rnodes){
  
  SEXP result = NILSXP;
  result = PROTECT(Rf_allocVector(VECSXP, 2));
  
  {
    //split1D -> 0: splitted_edges
    //        -> 1: vectors of global midpoints numbers
    split1D(result, Rnodes, Redges, 0);
    // midpoints coordinates
    compute_midpoints(result, Rnodes, Redges, 1);
  }
  
  UNPROTECT(1);
  
  return result;
  }

//! A function to refine a 1.5D mesh.
/*!
Each edge is splitted into subedges whose lengths are less or equal to the maximum allowed length \p Rdelta
\param Rnodes an R-matrix storing the mesh nodes
 \param Redges an R-matrix storing the mesh edges
 \param Rdelta an R integer representing the maximum allowed length
*/
SEXP refine1D(SEXP Rnodes, SEXP Redges, SEXP Rdelta){

    const RIntegerMatrix edges_old(Redges);
    const RNumericMatrix nodes(Rnodes);
    Real delta = REAL(Rdelta)[0];
    // stores the length of the i-th old edge;
    std::vector<Real> lengths(edges_old.nrows(),0);
    // stores the number of sub edges in which the i-th old edge is splitted into
    std::vector<UInt> num_subs(edges_old.nrows(),0);

    UInt num_new_nodes = 0;
    UInt num_tot_edges = 0;
    for(UInt i=0; i<edges_old.nrows(); ++i)
    {
        Real x0 = nodes(edges_old(i,0),0);
        Real y0 = nodes(edges_old(i,0),1);
        Real x1 = nodes(edges_old(i,1),0);
        Real y1 = nodes(edges_old(i,1),1);

        lengths[i] = std::sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );

        // number of sub_edges
        num_subs[i] = lengths[i] > delta ? std::ceil( lengths[i]/delta ) : 1;

        // number of new nodes
        num_new_nodes += num_subs[i] - 1;
        // number of edges
        num_tot_edges += num_subs[i];
    }

    SEXP result = PROTECT(Rf_allocVector(VECSXP,2));

    SET_VECTOR_ELT(result, 1, Rf_allocMatrix(INTSXP, num_tot_edges, 2));
    RIntegerMatrix edges( VECTOR_ELT(result,1));

    SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, num_new_nodes,2));
    RNumericMatrix nodes_new(VECTOR_ELT(result,0));

    UInt total_nodes = nodes.nrows();
    UInt edges_count = 0;
    UInt nodes_new_count = 0;

    static constexpr Real eps = 10 * std::numeric_limits<Real>::epsilon();

    for( UInt i=0; i<edges_old.nrows(); ++i)
    {
        if( num_subs[i] == 1)
        {
            // Indexes in R starts from 1, in C++ from 0, needed transformations!
            edges( edges_count,0) = edges_old(i,0) + 1;
            edges( edges_count,1) = edges_old(i,1) + 1;
            ++edges_count;
        }
        else // there is at least one new internal node!
        {
            Real x0 = nodes(edges_old(i,0),0);
            Real y0 = nodes(edges_old(i,0),1);
            Real x1 = nodes(edges_old(i,1),0);
            Real y1 = nodes(edges_old(i,1),1);

            Real cos_;
            Real sin_;

            if( std::abs(x1-x0) < eps && (y1-y0) > 0  ){
                cos_ =  0.;
                sin_ =  1.;
            }else if( std::abs( x1-x0) < eps && (y1-y0) < 0 ){
                cos_ =  0.;
                sin_ = -1.;
            }else if( std::abs( y1-y0)<eps && (x1-x0) > 0 ) {
                cos_ =  1.;
                sin_ =  0.;
            }else if( std::abs( y1-y0)<eps && (x1-x0) < 0 ){
                cos_ = -1.;
                sin_ =  0.;
            }else{
                Real slope = (y1-y0)/(x1-x0);
                cos_ = (x1-x0)>0. ? std::sqrt(1./(1.+slope*slope)) : -std::sqrt(1./(1.+slope*slope)) ;
                sin_ = (y1-y0)>0. ? std::sqrt(slope*slope/(1.+slope*slope)) : -std::sqrt(slope*slope/(1.+slope*slope)) ;
            }

            Real delta_ = lengths[i]/num_subs[i];
            // vector storing the global numbers of the new nodes which belong to old current edges
            // If N = num_subs[i] then there are N-1 internal nodes (the newest ones)
            std::vector<UInt> nodes_global_num( num_subs[i] + 1,0);

            //Fixing the first and the last nodes
            nodes_global_num[0] = edges_old(i,0);
            nodes_global_num[ num_subs[i] ] = edges_old(i,1);

            // Computing new nodes coordinates
            for(UInt n=0; n < num_subs[i]-1 ; ++n)
            {
                nodes_global_num[n + 1] = n + total_nodes;
                nodes_new(nodes_new_count + n,0) = x0 + delta_ * cos_ * (n + 1);
                nodes_new(nodes_new_count + n,1) = y0 + delta_ * sin_ * (n + 1);
            }

            nodes_new_count += num_subs[i] - 1;
            total_nodes+= num_subs[i] - 1;

            for(UInt e=0; e<num_subs[i]; ++e)
            {
                // Indexes in R starts from 1, in C++ from 0, needed transformations!
                edges(edges_count,0) =  nodes_global_num[e] + 1;
                edges(edges_count,1) = nodes_global_num[e+1] + 1;
                ++edges_count;
            }
        }

    }

    UNPROTECT(1);
    return result;
 }

}
