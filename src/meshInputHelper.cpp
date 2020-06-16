#define R_VERSION_

#include "fdaPDE.h"

#include "mesh_input_helper.h"

#include <array>

extern "C" {


SEXP CPP_SurfaceMeshHelper(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};
  using OutputType = typename simplex_container<2>::OutputType;

  int *triangles = INTEGER(Rtriangles);

  UInt num_triangles = INTEGER(Rf_getAttrib(Rtriangles, R_DimSymbol))[0];
  UInt num_nodes = INTEGER(Rf_getAttrib(Rnodes, R_DimSymbol))[0];

  OutputType out;

  {
  simplex_container<2> edges_list(triangles, num_triangles, num_nodes, EDGES_ORDERING);

  out=edges_list.assemble_output();
  }

  const std::vector<UInt> &edges=std::get<0>(out);
  const std::vector<bool> &edgesmarkers=std::get<1>(out);
  const std::vector<bool> &nodesmarkers=std::get<2>(out);
  const std::vector<int> &neighbors=std::get<3>(out);

  SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 4));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, edges.size()/2, 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(LGLSXP, edgesmarkers.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(LGLSXP, nodesmarkers.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocMatrix(INTSXP, neighbors.size()/3, 3));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
  for (UInt i=0; i<edges.size(); ++i)
    rans[i] = edges[i]+1;

	int *rans1 = LOGICAL(VECTOR_ELT(result, 1));
  for (UInt i=0; i<edgesmarkers.size(); ++i)
    rans1[i] = edgesmarkers[i];

	int *rans2 = LOGICAL(VECTOR_ELT(result, 2));
  for (UInt i=0; i<nodesmarkers.size(); ++i)
    rans2[i] = nodesmarkers[i];

	int *rans3 = INTEGER(VECTOR_ELT(result, 3));
  for (UInt i=0; i<neighbors.size(); ++i)
    rans3[i] = neighbors[i]+1;


	UNPROTECT(1);

  return result;
}


SEXP CPP_SurfaceMeshOrder2(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};
  using OutputType = typename simplex_container<2>::OutputType;

  int *triangles = INTEGER(Rtriangles);
  double *nodes = REAL(Rnodes);

  UInt num_triangles = INTEGER(Rf_getAttrib(Rtriangles, R_DimSymbol))[0];
  UInt num_nodes = INTEGER(Rf_getAttrib(Rnodes, R_DimSymbol))[0];

  OutputType out;
  std::vector<int> extended_triangles;

  {
  simplex_container<2> edges_list(triangles, num_triangles, num_nodes, EDGES_ORDERING);

  out=edges_list.assemble_output();
  extended_triangles=order2extend(edges_list);
  }

  const std::vector<UInt> &edges=std::get<0>(out);
  const std::vector<bool> &edgesmarkers=std::get<1>(out);
  std::vector<bool> &nodesmarkers=std::get<2>(out);
  const std::vector<int> &neighbors=std::get<3>(out);

  std::vector<double> midpoints{compute_midpoints(nodes, edges, num_nodes)};

  nodesmarkers.resize(nodesmarkers.size()+midpoints.size()/3, false);


  SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 6));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, edges.size()/2, 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(LGLSXP, edgesmarkers.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(LGLSXP, nodesmarkers.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocMatrix(INTSXP, neighbors.size()/3, 3));
  SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, midpoints.size()/3, 3));
  SET_VECTOR_ELT(result, 5, Rf_allocMatrix(INTSXP, extended_triangles.size()/3, 3));


	int *rans = INTEGER(VECTOR_ELT(result, 0));
  for (UInt i=0; i<edges.size(); ++i)
    rans[i] = edges[i]+1;

	int *rans1 = LOGICAL(VECTOR_ELT(result, 1));
  for (UInt i=0; i<edgesmarkers.size(); ++i)
    rans1[i] = edgesmarkers[i];

	int *rans2 = LOGICAL(VECTOR_ELT(result, 2));
  for (UInt i=0; i<nodesmarkers.size(); ++i)
    rans2[i] = nodesmarkers[i];

	int *rans3 = INTEGER(VECTOR_ELT(result, 3));
  for (UInt i=0; i<neighbors.size(); ++i)
    rans3[i] = neighbors[i]+1;

  double *rans4 = REAL(VECTOR_ELT(result, 4));
  for (UInt i=0; i<midpoints.size(); ++i)
    rans4[i] = midpoints[i];

  int *rans5 = INTEGER(VECTOR_ELT(result, 5));
  for (UInt i=0; i<extended_triangles.size(); ++i)
    rans5[i] = extended_triangles[i];


	UNPROTECT(1);

  return result;
}

SEXP CPP_TriangleMeshSplit(SEXP Rtriangles, SEXP Rnodes){

  static constexpr std::array<UInt, 6> EDGES_ORDERING = {1,2,0,2,0,1};

  int *triangles = INTEGER(Rtriangles);
  double *nodes = REAL(Rnodes);


  UInt num_triangles = INTEGER(Rf_getAttrib(Rtriangles, R_DimSymbol))[0];
  UInt num_nodes = INTEGER(Rf_getAttrib(Rnodes, R_DimSymbol))[0];

  std::vector<double> midpoints;
  std::vector<UInt> splitted_triangles;

  {
  simplex_container<2> edges_list(triangles, num_triangles, num_nodes, EDGES_ORDERING);

  std::vector<UInt> extended_triangles{order2extend(edges_list)};
  std::vector<UInt> edges{edges_list.get_simplexes()};

  midpoints=compute_midpoints(nodes, edges, num_nodes);
  splitted_triangles=split(extended_triangles, triangles, num_triangles);
  }

  SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, splitted_triangles.size()/3, 3));
  SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, midpoints.size()/3, 3));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
  for (UInt i=0; i<splitted_triangles.size(); ++i)
    rans[i] = splitted_triangles[i];

  double *rans1 = REAL(VECTOR_ELT(result, 1));
  for (UInt i=0; i<midpoints.size(); ++i)
    rans1[i] = midpoints[i];


	UNPROTECT(1);

  return result;
}


SEXP CPP_VolumeMeshHelper(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> FACES_ORDERING = {1,2,3,0,2,3,0,1,3,0,1,2};
  using OutputType = typename simplex_container<3>::OutputType;


  int *tetrahedrons = INTEGER(Rtetrahedrons);

  UInt num_tetrahedrons = INTEGER(Rf_getAttrib(Rtetrahedrons, R_DimSymbol))[0];
  UInt num_nodes = INTEGER(Rf_getAttrib(Rnodes, R_DimSymbol))[0];

  OutputType out;

  {
  simplex_container<3> faces_list(tetrahedrons, num_tetrahedrons, num_nodes, FACES_ORDERING);

  out=faces_list.assemble_output();
  }

  const std::vector<UInt> &faces=std::get<0>(out);
  const std::vector<bool> &facesmarkers=std::get<1>(out);
  const std::vector<bool> &nodesmarkers=std::get<2>(out);
  const std::vector<int> &neighbors=std::get<3>(out);

  SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 4));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, faces.size()/3, 3));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(LGLSXP, facesmarkers.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(LGLSXP, nodesmarkers.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocMatrix(INTSXP, neighbors.size()/4, 4));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
  for (UInt i=0; i<faces.size(); ++i)
    rans[i] = faces[i]+1;

	int *rans1 = LOGICAL(VECTOR_ELT(result, 1));
  for (UInt i=0; i<facesmarkers.size(); ++i)
    rans1[i] = facesmarkers[i];

	int *rans2 = LOGICAL(VECTOR_ELT(result, 2));
  for (UInt i=0; i<nodesmarkers.size(); ++i)
    rans2[i] = nodesmarkers[i];

	int *rans3 = INTEGER(VECTOR_ELT(result, 3));
  for (UInt i=0; i<neighbors.size(); ++i)
    rans3[i] = neighbors[i]+1;


	UNPROTECT(1);

  return result;

}

SEXP CPP_VolumeMeshOrder2(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> FACES_ORDERING = {1,2,3,0,2,3,0,1,3,0,1,2};
  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};
  using OutputType = typename simplex_container<3>::OutputType;

  int *tetrahedrons = INTEGER(Rtetrahedrons);
  double *nodes = REAL(Rnodes);

  UInt num_tetrahedrons = INTEGER(Rf_getAttrib(Rtetrahedrons, R_DimSymbol))[0];
  UInt num_nodes = INTEGER(Rf_getAttrib(Rnodes, R_DimSymbol))[0];

  OutputType out;

  {
  simplex_container<3> faces_list(tetrahedrons, num_tetrahedrons, num_nodes, FACES_ORDERING);

  out=faces_list.assemble_output();
  }

  const std::vector<UInt> &faces=std::get<0>(out);
  const std::vector<bool> &facesmarkers=std::get<1>(out);
  std::vector<bool> &nodesmarkers=std::get<2>(out);
  const std::vector<int> &neighbors=std::get<3>(out);


  std::vector<int> extended_tetrahedrons;
  std::vector<double> midpoints;

  {
  simplex_container<2> edges_list(tetrahedrons, num_tetrahedrons, num_nodes, EDGES_ORDERING);

  extended_tetrahedrons=order2extend(edges_list);
  std::vector<UInt> edges{edges_list.get_simplexes()};
  midpoints=compute_midpoints(nodes, edges, num_nodes);
  }

  nodesmarkers.resize(nodesmarkers.size()+midpoints.size()/3, false);


  SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 4));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, faces.size()/3, 3));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(LGLSXP, facesmarkers.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(LGLSXP, nodesmarkers.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocMatrix(INTSXP, neighbors.size()/4, 4));
  SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, midpoints.size()/3, 3));
  SET_VECTOR_ELT(result, 5, Rf_allocMatrix(INTSXP, extended_tetrahedrons.size()/6, 6));


	int *rans = INTEGER(VECTOR_ELT(result, 0));
  for (UInt i=0; i<faces.size(); ++i)
    rans[i] = faces[i]+1;

	int *rans1 = LOGICAL(VECTOR_ELT(result, 1));
  for (UInt i=0; i<facesmarkers.size(); ++i)
    rans1[i] = facesmarkers[i];

	int *rans2 = LOGICAL(VECTOR_ELT(result, 2));
  for (UInt i=0; i<nodesmarkers.size(); ++i)
    rans2[i] = nodesmarkers[i];

	int *rans3 = INTEGER(VECTOR_ELT(result, 3));
  for (UInt i=0; i<neighbors.size(); ++i)
    rans3[i] = neighbors[i]+1;

  double *rans4 = REAL(VECTOR_ELT(result, 4));
  for (UInt i=0; i<midpoints.size(); ++i)
    rans4[i] = midpoints[i];

  int *rans5 = INTEGER(VECTOR_ELT(result, 5));
  for (UInt i=0; i<extended_tetrahedrons.size(); ++i)
    rans5[i] = extended_tetrahedrons[i];


	UNPROTECT(1);

  return result;

}


SEXP CPP_TetraMeshSplit(SEXP Rtetrahedrons, SEXP Rnodes){

  static constexpr std::array<UInt, 12> EDGES_ORDERING = {0,1,0,2,0,3,1,2,2,3,1,3};

  int *tetrahedrons = INTEGER(Rtetrahedrons);
  double *nodes = REAL(Rnodes);

  UInt num_tetrahedrons = INTEGER(Rf_getAttrib(Rtetrahedrons, R_DimSymbol))[0];
  UInt num_nodes = INTEGER(Rf_getAttrib(Rnodes, R_DimSymbol))[0];

  std::vector<double> midpoints;
  std::vector<UInt> splitted_tetrahedrons;

  {
  simplex_container<2> edges_list(tetrahedrons, num_tetrahedrons, num_nodes, EDGES_ORDERING);

  std::vector<int> extended_tetrahedrons{order2extend(edges_list)};
  std::vector<UInt> edges{edges_list.get_simplexes()};

  midpoints=compute_midpoints(nodes, edges, num_nodes);
  splitted_tetrahedrons=split3D(extended_tetrahedrons, tetrahedrons, num_tetrahedrons);
  }


  SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, splitted_tetrahedrons.size()/4, 4));
  SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, midpoints.size()/3, 3));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
  for (UInt i=0; i<splitted_tetrahedrons.size(); ++i)
    rans[i] = splitted_tetrahedrons[i];

  double *rans1 = REAL(VECTOR_ELT(result, 1));
  for (UInt i=0; i<midpoints.size(); ++i)
    rans1[i] = midpoints[i];


	UNPROTECT(1);

  return result;
}


}
