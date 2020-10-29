#ifndef __POINT_H__
#define __POINT_H__

#include "../../FdaPDE.h"

/*! \brief A class storing an ID and providing the relevant methods
 *
 *  This class is inherited from Point and Element classes
 *  It is used to store a UInt representing an ID and it provides
 *  some methods to check if the ID is valid
 */

struct Identifier{

  //! The biggest number of type UInt represents a not valid ID (NVAL)
	static constexpr UInt NVAL = std::numeric_limits<UInt>::max();

  //! A method to check if the object ID has been assigned
	bool unassignedId() const {return id_==NVAL;}
  //! A method to check if the object ID has been assigned
  bool hasValidId() const {return !unassignedId();}
  //! A method to check if the object BCId has been assigned
  /*
   * This method is no longer relevant in the present version of the code and could be removed
   */
	bool unassignedBc() const {return bcId_==NVAL;}

  //! A method to get the object ID
	const UInt& id() const {return id_;}
  //! A method to get the object BCId
  /*
   * This method is no longer relevant in the present version of the code and could be removed
   */
	const UInt& bcId() const {return bcId_;}
  //! A method to get the object ID
	const UInt& getId() const {return id_;}

protected:
  //! A UInt storing the object ID
  /*
   * The default object ID is not valid if not explicitly set
   */
  UInt id_=NVAL;
  // Note: bcId was probably intended to mark boundary objects but is never actually
  // used, so it could be removed
  //! A UInt storing the object BCId
  /*
   * It was probably intended to mark boundary objects but is never actually used and it could be removed
   */
	UInt bcId_=NVAL;

  //! The default constructor
  /*
   * The default constructor initializes ID as not valid
   * It is constexpr so that the class can be used as a literal type
   * It is protected to make the class act as an abstract base class
   * without any runtime cost due to virtual functions
   */
  constexpr Identifier()=default;

  //! A constructor explicitly setting the ID
  /*
   * This constructor initializes ID using the given UInt value
   * It is constexpr so that the class can be used as a literal type
   * It is protected to make the class act as an abstract base class
   * without any runtime cost due to virtual functions
   * \param id a UInt representing the ID to be set
   */
  constexpr Identifier(UInt id) : id_(id) {}

  //! A constructor explicitly setting ID and BCId
  /*
   * This constructor initializes ID and BCId using the given UInt value
   * Since BCId is never actually used this constructor could be removed
   * It is constexpr so that the class can be used as a literal type
   * It is protected to make the class act as an abstract base class
   * without any runtime cost due to virtual functions
   * \param id a UInt representing the ID to be set
   * \param bcId a UInt representing the BCId to be set
   */
  constexpr Identifier(UInt id, UInt bcId) : id_(id), bcId_(bcId) {}

  //!The default destructor
  /*
   * It is protected to be consistent with the constructor being protected, so that
   * it is not possible to delete an object inheriting from Identifier through a pointer to Identifier
   */
  ~Identifier()=default;
};



template<UInt ndim>
class Point : public Identifier{
  static_assert(ndim==2 || ndim==3,
								 "ERROR! TRYING TO INSTANTIATE POINT IN UNIMPLEMENTED DIMENSION! See mesh_objects.h");
  public:
    using pointCoords = std::array<Real,ndim>;
    using EigenCoords = Eigen::Matrix<Real,ndim,1>;
    using EigenMap2Coords = Eigen::Map<EigenCoords>;
    using EigenMap2ConstCoords = Eigen::Map<const EigenCoords>;


		// Note: these don't really mean anything, they're just here for compatibility
		// with the adtree implementation
		static constexpr UInt dp() {return 3;}
		static constexpr UInt dt() {return 3;}
		static constexpr UInt coordsize() {return 3;}

    // Note: some of the constructors are declared constexpr so that Point can be
    // used as a literal type (see integration.h)
    //! The default constructor initializing the origin
    constexpr Point() : coord_() {}
    // Full constructor initializing every member explicitly
    constexpr Point(UInt id, UInt bcId, const pointCoords& coord) :
              Identifier(id, bcId), coord_(coord) {}
    // Additional constructors initializing only some members
    constexpr Point(const pointCoords& coord) :
                Point(NVAL, NVAL, coord) {}
    constexpr Point(UInt id, const pointCoords& coord) :
                Point(id, NVAL, coord) {}
    // Additional constructor (needed for disambiguation when constructing from
    // bracketed initializer lists)
    constexpr Point(UInt id, UInt bcId, const Real(&coord)[ndim]);
    constexpr Point(const Real(&coord)[ndim]) :
                Point(NVAL, NVAL, coord) {}
    constexpr Point(UInt id, const Real(&coord)[ndim]) :
                Point(id, NVAL, coord) {}

    // Additional constructor allowing construction from an eigen vector
    Point(UInt id, UInt bcId, const EigenCoords &coord);
    Point(const EigenCoords &coord) :
      Point(NVAL, NVAL, coord) {}
    Point(UInt id, const EigenCoords &coord) :
      Point(id, NVAL, coord) {}

    // Additional constructor for convenience in dealing with R data (e.g. meshes)
    Point(UInt id, const RNumericMatrix& points);

    // Views!
    EigenMap2Coords eigenView() {return EigenMap2Coords(&coord_[0]);}
    EigenMap2ConstCoords eigenView() const {return EigenMap2ConstCoords(&coord_[0]);}

    EigenMap2ConstCoords eigenConstView() const {return EigenMap2ConstCoords(&coord_[0]);}


    // Overloaded subscript operator
    Real& operator[](UInt i) {return coord_[i];}
    const Real& operator[](UInt i) const {return coord_[i];}

		pointCoords& coord() {return coord_;}
		const pointCoords& coord() const {return coord_;}

    // Member functions returning the distance/squared distance between two points
    // Note: both are defined because the squared distance function is more efficient
    // and should be preferred if one is only interested in comparing distances
    Real dist(const Point&) const;
    Real dist2(const Point&) const;
    // Friend functions returning the distance/squared distance between two points
    friend Real dist(const Point& lhs, const Point& rhs) {return lhs.dist(rhs);};
    friend Real dist2(const Point& lhs, const Point& rhs) {return lhs.dist2(rhs);};

    // Overloaded "+="/"-=" operators
    // These operators add/subtract the coordinates elementwise
    Point& operator+=(const Point&);
    Point& operator-=(const Point&);
    // Overloaded "+"/"-" operator
    // These operators add/subtract the coordinates elementwise and return a new point
    friend Point operator+(Point lhs, const Point& rhs) {return lhs+=rhs;};
    friend Point operator-(Point lhs, const Point& rhs) {return lhs-=rhs;};

    //! Overload the << operator to easily print Point info (note: define it in class
    // to avoid a forward declaration)
    template <UInt NDIM>
    friend std::ostream& operator<<(std::ostream& os, const Point<NDIM>& p);

  private:

    pointCoords coord_;

  };


#include "Point_imp.h"

#endif
