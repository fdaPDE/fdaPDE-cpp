#ifndef ROBJECTS_H
#define ROBJECTS_H


class RNumericMatrix{
public:

	RNumericMatrix(Real * const matr, const UInt nrows, const UInt ncols) :
		matr_(matr), nrows_(nrows), ncols_(ncols) {}

	RNumericMatrix(SEXP matr) : 
		matr_(REAL(matr)), 
			nrows_(INTEGER(Rf_getAttrib(matr, R_DimSymbol))[0]), 
				ncols_(INTEGER(Rf_getAttrib(matr, R_DimSymbol))[1]) {}

	Real& operator()(UInt i, UInt j) {return matr_[i+nrows_*j];}
	const Real& operator()(UInt i, UInt j) const {return matr_[i+nrows_*j];}
	
	Real& operator[](UInt i) {return matr_[i];}
	const Real& operator[](UInt i) const {return matr_[i];}

	const UInt& nrows() const {return nrows_;}
	const UInt& ncols() const {return ncols_;}

private:
	Real * const matr_;
	const UInt nrows_;
	const UInt ncols_;

};

class RIntegerMatrix{
public:

	RIntegerMatrix(UInt * const matr, const UInt nrows, const UInt ncols) :
		matr_(matr), nrows_(nrows), ncols_(ncols) {}


	RIntegerMatrix(SEXP matr) : 
		matr_(INTEGER(matr)), 
			nrows_(INTEGER(Rf_getAttrib(matr, R_DimSymbol))[0]), 
				ncols_(INTEGER(Rf_getAttrib(matr, R_DimSymbol))[1]) {}
	

	UInt& operator()(UInt i, UInt j) {return matr_[i+nrows_*j];}
	const UInt& operator()(UInt i, UInt j) const {return matr_[i+nrows_*j];}
	
	UInt& operator[](UInt i) {return matr_[i];}
	const UInt& operator[](UInt i) const {return matr_[i];}

	const UInt& nrows() const {return nrows_;}
	const UInt& ncols() const {return ncols_;}

private:
	UInt * const matr_;
	const UInt nrows_;
	const UInt ncols_;

};

// A class to access/store matrix of matrix
// It is used to store neighbors info in 1D mesh
class RIntMatrixMatrix{
public:
    using value_type = RIntegerMatrix;
    using container_type = std::vector< value_type >;

    RIntMatrixMatrix(SEXP Robject):nrows_( INTEGER(Rf_getAttrib(Robject, R_DimSymbol))[0] ),
                                   ncols_( INTEGER(Rf_getAttrib(Robject, R_DimSymbol))[1] ){
        matr_.reserve(nrows_*ncols_);
        for(UInt i=0; i<nrows_*ncols_; ++i){
            matr_.emplace_back(VECTOR_ELT(Robject,i));
        }
    }

    value_type& operator()(UInt i , UInt j){  return matr_[i+nrows_*j];}
    const value_type& operator() (UInt i, UInt j) const{ return matr_[i+nrows_*j];}

    value_type& operator[](UInt j){ return matr_[j];}
    const value_type& operator[] (UInt j) const{return matr_[j];}

    const UInt& nrows() const { return nrows_;}
    const UInt& ncols() const { return ncols_;}

private:
    container_type matr_;
    const UInt nrows_;
    const UInt ncols_;

};

#endif
