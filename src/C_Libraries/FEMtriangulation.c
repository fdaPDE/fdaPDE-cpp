#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdint.h>
#undef PI
#define printf Rprintf
//#define VOID void
#define TRILIBRARY
//#define ANSI_DECLARATORS

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define TRIREAL float
#else /* not SINGLE */
#define TRIREAL double
#endif /* not SINGLE */



#include "triangle.h"

/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/

//void report(io, markers, reporttriangles, reportneighbors, reportsegments,
//            reportedges, reportnorms)
//{
//struct triangulateio *io;
//int markers;
//int reporttriangles;
//int reportneighbors;
//int reportsegments;
//int reportedges;
//int reportnorms;
//
//  int i, j;
//
//  for (i = 0; i < io->numberofpoints; i++) {
//    printf("Point %4d:", i);
//    for (j = 0; j < 2; j++) {
//      printf("  %.6g", io->pointlist[i * 2 + j]);
//    }
//    if (io->numberofpointattributes > 0) {
//      printf("   attributes");
//    }
//    for (j = 0; j < io->numberofpointattributes; j++) {
//      printf("  %.6g",
//             io->pointattributelist[i * io->numberofpointattributes + j]);
//    }
//    if (markers) {
//      printf("   marker %d\n", io->pointmarkerlist[i]);
//    } else {
//      printf("\n");
//    }
//  }
//  printf("\n");
//
//  if (reporttriangles || reportneighbors) {
//    for (i = 0; i < io->numberoftriangles; i++) {
//      if (reporttriangles) {
//        printf("Triangle %4d points:", i);
//        for (j = 0; j < io->numberofcorners; j++) {
//          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
//        }
//        if (io->numberoftriangleattributes > 0) {
//          printf("   attributes");
//        }
//        for (j = 0; j < io->numberoftriangleattributes; j++) {
//          printf("  %.6g", io->triangleattributelist[i *
//                                         io->numberoftriangleattributes + j]);
//        }
//        printf("\n");
//      }
//      if (reportneighbors) {
//        printf("Triangle %4d neighbors:", i);
//        for (j = 0; j < 3; j++) {
//          printf("  %4d", io->neighborlist[i * 3 + j]);
//        }
//        printf("\n");
//      }
//    }
//    printf("\n");
//  }
//
//  if (reportsegments) {
//    for (i = 0; i < io->numberofsegments; i++) {
//      printf("Segment %4d points:", i);
//      for (j = 0; j < 2; j++) {
//        printf("  %4d", io->segmentlist[i * 2 + j]);
//      }
//      if (markers) {
//        printf("   marker %d\n", io->segmentmarkerlist[i]);
//      } else {
//        printf("\n");
//      }
//    }
//    printf("\n");
//  }
//
//  if (reportedges) {
//    for (i = 0; i < io->numberofedges; i++) {
//      printf("Edge %4d points:", i);
//      for (j = 0; j < 2; j++) {
//        printf("  %4d", io->edgelist[i * 2 + j]);
//      }
//      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
//        for (j = 0; j < 2; j++) {
//          printf("  %.6g", io->normlist[i * 2 + j]);
//        }
//      }
//      if (markers) {
//        printf("   marker %d\n", io->edgemarkerlist[i]);
//      } else {
//        printf("\n");
//      }
//    }
//    printf("\n");
//  }
//}

/*****************************************************************************/
/*                                                                           */
/*  R_triangulate_native()   Function that allows a "native" usage of the library            */
/*                                                                           */
/*****************************************************************************/

SEXP R_triangulate_native(SEXP P, SEXP PB, SEXP PA, SEXP S, SEXP SB, SEXP(H), SEXP T, SEXP Rflags)
{
  /* Output variables */
  SEXP oP, oPB, oPA, oT, oS, oSB, oE, oEB, oNT, oPV, oEV, oNV, oAV;
  SEXP ans;
  double *xoP, *xoPA, *xoPV, *xoNV, *xoAV;
  int *xoT, *xoPB, *xoS, *xoSB, *xoE, *xoEB, *xoEV, *xoNT;

  const char *flags = CHAR(STRING_ELT(Rflags,0));

  //printf("Flag string is: %s \n", flags);

  /* Convert input point matrix into array */
  PROTECT(P = AS_NUMERIC(P));
  /* Convert input boundary markers into array */
  // PROTECT(B = AS_NUMERIC(B));

  /* Create the triangulateio structures */
  struct triangulateio in, mid, vorout;

  in.numberofpoints = LENGTH(P)/2;
  in.pointlist = REAL(P);

  if(Rf_length(PB) == 0)
	  in.pointmarkerlist = NULL;
  else
	  in.pointmarkerlist = INTEGER(PB);

  in.numberofpointattributes = Rf_ncols(PA);
  in.pointattributelist = REAL(PA);

  in.numberofsegments = LENGTH(S)/2;
  in.segmentlist = INTEGER(S);
  if(Rf_length(SB) == 0)
  	in.segmentmarkerlist = NULL;
  else
    in.segmentmarkerlist = INTEGER(SB);
  in.numberofholes = LENGTH(H)/2;
  in.holelist = REAL(H);				/* Not needed if -E switch used. */

  in.numberofregions = 0;

  in.numberoftriangles 			= INTEGER(getAttrib(T, R_DimSymbol))[1];
  in.numberofcorners 			= INTEGER(getAttrib(T, R_DimSymbol))[0];
  in.trianglelist 				= INTEGER(T);          /* Not needed if -E switch used. */

//  in.trianglelist = (int *) NULL;
  in.triangleattributelist 		= (TRIREAL *) NULL;  //we do not consider this option
  in.trianglearealist 			= (TRIREAL *) NULL;
  //in.numberoftriangles 			= 0;
  //in.numberofcorners 			= 3;
  in.numberoftriangleattributes = 0;


  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */

  mid.pointlist = (TRIREAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  mid.pointattributelist = (TRIREAL *) NULL;
  mid.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  mid.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  mid.triangleattributelist = (TRIREAL *) NULL;
  mid.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  /* Needed only if segments are output (-p or -c) and -P not used: */
  mid.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  mid.segmentmarkerlist = (int *) NULL;
  mid.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  mid.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
  mid.holelist 		 = (TRIREAL *) NULL;



  vorout.pointlist = (TRIREAL *) NULL;        /* Needed only if -v switch used. */
  /* Needed only if -v switch used and number of attributes is not zero: */
  vorout.pointattributelist = (TRIREAL *) NULL;
  vorout.edgelist = (int *) NULL;          /* Needed only if -v switch used. */
  vorout.normlist = (TRIREAL *) NULL;         /* Needed only if -v switch used. */

  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), preserve the convex hull (c), assign a regional       */
  /*   attribute to each element (A), and                              */
  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  /*   neighbor list (n).                                              */

  //char flags[200];
  //strcpy(flags, inputflags);


  triangulate(strdup(flags), &in, &mid, &vorout);

//printf("Initial triangulation:\n\n");
      //report(&mid, 1, 1, 1, 1, 1, 0);
  /* printf("Initial Voronoi diagram:\n\n");
     report(&vorout, 0, 0, 0, 0, 1, 1); */

  /* Attach area constraints to the triangles in preparation for */
  /*   refining the triangulation.                               */

  /* /\* Needed only if -r and -a switches used: *\/ */
   //mid.trianglearealist = (TRIREAL *) malloc(mid.numberoftriangles * sizeof(TRIREAL));
   //mid.trianglearealist[0] = 3.0;
   //mid.trianglearealist[1] = 1.0;

  /* /\* Make necessary initializations so that Triangle can return a *\/ */
  /* /\*   triangulation in `out'.                                    *\/ */

  /* out.pointlist = (TRIREAL *) NULL;            /\* Not needed if -N switch used. *\/ */
  /* /\* Not needed if -N switch used or number of attributes is zero: *\/ */
  /* out.pointattributelist = (TRIREAL *) NULL; */
  /* out.trianglelist = (int *) NULL;          /\* Not needed if -E switch used. *\/ */
  /* /\* Not needed if -E switch used or number of triangle attributes is zero: *\/ */
  /* out.triangleattributelist = (TRIREAL *) NULL; */

  /* /\* Refine the triangulation according to the attached *\/ */
  /* /\*   triangle area constraints.                       *\/ */

  /* triangulate("praBP", &mid, &out, (struct triangulateio *) NULL); */

  /* printf("Refined triangulation:\n\n"); */
  /* report(&out, 0, 1, 0, 0, 0, 0); */

  /* Make space for answers */
  PROTECT(oP  = allocMatrix(REALSXP,  mid.numberofpoints, 2));
  PROTECT(oPB = allocMatrix(INTSXP,   mid.numberofpoints, 1));
  PROTECT(oPA = allocMatrix(REALSXP,  mid.numberofpoints, mid.numberofpointattributes));
  PROTECT(oT  = allocMatrix(INTSXP,   mid.numberoftriangles, mid.numberofcorners));
  PROTECT(oS  = allocMatrix(INTSXP,   mid.numberofsegments, 2));
  PROTECT(oSB = allocMatrix(INTSXP,   mid.numberofsegments, 1));
  PROTECT(oE  = allocMatrix(INTSXP,   mid.numberofedges, 2));
  PROTECT(oEB = allocMatrix(INTSXP,   mid.numberofedges, 1));
  PROTECT(oNT = allocMatrix(INTSXP,   mid.numberoftriangles, 3));
  PROTECT(oPV  = allocMatrix(REALSXP, vorout.numberofpoints, 2));
  PROTECT(oEV  = allocMatrix(INTSXP,  vorout.numberofedges, 2));
  PROTECT(oNV  = allocMatrix(REALSXP, vorout.numberofpoints, 2));
  PROTECT(oAV  = allocMatrix(REALSXP, vorout.numberofpoints, mid.numberofpointattributes));

  //printf("The number of corners is %d", mid.numberofcorners);

  xoP = REAL(oP);
  for (int i = 0; i < mid.numberofpoints; i++) {
    for (int j = 0; j < 2; j++) {
      xoP[j * mid.numberofpoints + i] = mid.pointlist[i * 2 + j];
    }
  }

  xoPB = INTEGER(oPB);
  for (int i = 0; i < mid.numberofpoints; i++) {
    xoPB[i] = mid.pointmarkerlist[i];
  }

  xoPA = REAL(oPA);
  for (int i = 0; i < mid.numberofpoints; i++) {
    for (int j = 0; j < mid.numberofpointattributes; j++) {
      xoPA[j * mid.numberofpoints + i] = mid.pointattributelist[i * mid.numberofpointattributes + j];
    }
  }

  xoT = INTEGER(oT);
  for (int i = 0; i < mid.numberoftriangles; i++) {
    for (int j = 0; j < mid.numberofcorners; j++) {
      xoT[j * mid.numberoftriangles + i] = mid.trianglelist[i * mid.numberofcorners + j];
    }
  }

  xoS = INTEGER(oS);
  for (int i = 0; i < mid.numberofsegments; i++) {
    for (int j = 0; j < 2; j++) {
      xoS[j * mid.numberofsegments + i] = mid.segmentlist[i * 2 + j];
    }
  }

  xoSB = INTEGER(oSB);
  for (int i = 0; i < mid.numberofsegments; i++) {
    xoSB[i] = mid.segmentmarkerlist[i];
  }

  xoE = INTEGER(oE);
  for (int i = 0; i < mid.numberofedges; i++) {
    for (int j = 0; j < 2; j++) {
      xoE[j * mid.numberofedges + i] = mid.edgelist[i * 2 + j];
    }
  }

  xoEB = INTEGER(oEB);
  for (int i = 0; i < mid.numberofedges; i++) {
    xoEB[i] = mid.edgemarkerlist[i];
  }

  xoNT = INTEGER(oNT);
  for (int i = 0; i < mid.numberoftriangles; i++) {
    for (int j = 0; j < 3; j++) {
      xoNT[j * mid.numberoftriangles + i] = mid.neighborlist[i * 3 + j];
    }
  }

  xoPV = REAL(oPV);
  for (int i = 0; i < vorout.numberofpoints; i++) {
    for (int j = 0; j < 2; j++) {
      xoPV[j * vorout.numberofpoints + i] = vorout.pointlist[i * 2 + j];
    }
  }

  xoEV = INTEGER(oEV);
  for (int i = 0; i < vorout.numberofedges; i++) {
    for (int j = 0; j < 2; j++) {
      xoEV[j * vorout.numberofedges + i] = vorout.edgelist[i * 2 + j];
    }
  }

  xoNV = REAL(oNV);
  for (int i = 0; i < vorout.numberofpoints; i++) {
    for (int j = 0; j < 2; j++) {
      xoNV[j * vorout.numberofpoints + i] = vorout.normlist[i * 2 + j];
    }
  }

  xoAV = REAL(oAV);
  for (int i = 0; i < vorout.numberofpoints; i++) {
    for (int j = 0; j < mid.numberofpointattributes; j++) {
      xoAV[j * vorout.numberofpoints + i] = vorout.pointattributelist[i * mid.numberofpointattributes + j];
    }
  }

  PROTECT(ans = allocVector(VECSXP, 13));
  SET_VECTOR_ELT(ans, 0, oP);
  SET_VECTOR_ELT(ans, 1, oPB);
  SET_VECTOR_ELT(ans, 2, oPA);
  SET_VECTOR_ELT(ans, 3, oT);
  SET_VECTOR_ELT(ans, 4, oS);
  SET_VECTOR_ELT(ans, 5, oSB);
  SET_VECTOR_ELT(ans, 6, oE);
  SET_VECTOR_ELT(ans, 7, oEB);
  SET_VECTOR_ELT(ans, 8, oNT);
  SET_VECTOR_ELT(ans, 9, oPV);
  SET_VECTOR_ELT(ans, 10, oEV);
  SET_VECTOR_ELT(ans, 11, oNV);
  SET_VECTOR_ELT(ans, 12, oAV);
  UNPROTECT(15);

  /* Free all allocated arrays, including those allocated by Triangle. */
  Free(mid.pointlist);
  Free(mid.pointattributelist);
  Free(mid.pointmarkerlist);
  Free(mid.trianglelist);
  Free(mid.triangleattributelist);
  /* Free(mid.trianglearealist); */
  Free(mid.neighborlist);
  Free(mid.segmentlist);
  Free(mid.segmentmarkerlist);
  Free(mid.edgelist);
  Free(mid.edgemarkerlist);
  Free(vorout.pointlist);
  Free(vorout.pointattributelist);
  Free(vorout.edgelist);
  Free(vorout.normlist);
  /*  Free(out.pointlist);
  Free(out.pointattributelist);
  Free(out.trianglelist);
  Free(out.triangleattributelist); */

  return(ans);
}
