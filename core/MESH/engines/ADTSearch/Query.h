#ifndef __QUERY_H__
#define __QUERY_H__

template <unsigned int N> using rectangle = std::pair<SVector<N>, SVector<N>>;

// in geometric search a query is an N-dimensional rectangle described by its left-lower corner
// and upper-right corner. This class wraps a rectangle providing some facilities to perform geometric
// operations between queries
template <unsigned int N> class Query {
private:
  rectangle<N> queryRange_;
public:
  // constructor
  Query(const rectangle<N>& queryRange) : queryRange_(queryRange) {}

  // returns true if a given point lies inside the query
  bool contains(const SVector<N>& point) const;
  // returns true if the query intersects a given rectangle along a given dimension
  bool intersect(const rectangle<N>& rectangle) const; 

  rectangle<N> getRange() const { return queryRange_; }
};

// this function is called frequently during search
template <unsigned int N>
bool Query<N>::intersect(const rectangle<N> &rect) const {
  // keep track along which dimension query and rect intersects
  std::array<bool, N> boolVector{};
  
  for(size_t dim = 0; dim < N; ++dim){
    // load here once, access it fast
    double qs = queryRange_.second[dim], qf = queryRange_.first[dim];
    double rs = rect.second[dim],        rf = rect.first[dim];
    
    if(// partially overlapping sides
       (qs > rf && rs > qf) || (rs > qf && qs > rf) ||
       // queryRange_ contained into rect or viceversa 
       (rf < qf && qs < rs) || (qf < rf && rs < rs)){

      boolVector[dim] = true; // query intersect rectangle along this dimension
    }
  }

  // query and rect intersects if and only if they intersects along each dimension
  bool result = true;
  for(bool b : boolVector) result &= b;
  
  return result;
}

template <unsigned int N>
bool Query<N>::contains(const SVector<N> &point) const {
  // cycle over dimensions
  for(size_t dim = 0; dim < N; ++dim){
    // return false if point lies outside 
    if(point[dim] > queryRange_.second[dim] || point[dim] < queryRange_.first[dim])
      return false;
  }
  return true;
}

#endif // __QUERY_H__
