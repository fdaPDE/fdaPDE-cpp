#ifndef __MESH_INPUT_HELPER_IMP_HPP__
#define __MESH_INPUT_HELPER_IMP_HPP__


template<UInt mydim>
template<std::size_t SIZE>
void simplex_container<mydim>::fill_container(const UInt* const elements, const std::array<UInt, SIZE>& ORDERING){
 static_assert(SIZE==mydim*(mydim+1) || (mydim==2 && SIZE==12),
        "ERROR! ORDERING SIZE SHOULD BE EQUAL TO 2X THE NUMBER OF EDGES OR 3X THE NUMBER OF FACES! See: mesh_input_helper_imp.h");

 simplexes.reserve(num_elements*ORDERING.size()/mydim);

 {
   std::array<UInt,mydim> curr;
   for(UInt i=0; i<num_elements; ++i){
     for(UInt j=0; j<ORDERING.size()/mydim; ++j){
       for(UInt k=0; k<mydim; ++k)
        curr[k]=elements[i+num_elements*ORDERING[mydim*j+k]];
       std::sort(curr.begin(), curr.end());
       simplexes.emplace_back(simplex_t(i,j,curr));
     }
   }
 }

  bin_sort();
  check_duplicates();
  store_indexes();

}

template<UInt mydim>
void simplex_container<mydim>::bin_sort(){

  std::vector<UInt> positions;
  positions.reserve(simplexes.size());
  for(UInt i=0; i<simplexes.size(); ++i)
    positions.push_back(i);

  bin_sort_(mydim-1, positions);

  for(UInt i=0; i<positions.size(); ++i){
    UInt curr=i;
    while(i!=positions[curr]){
      UInt next=positions[curr];
      std::swap(simplexes[curr],simplexes[next]);
      positions[curr]=curr;
      curr=next;
    }
    positions[curr]=curr;
  }
}

// Recursive unction to sort container by ascending #(index+1) element of the arrays
template<UInt mydim>
void simplex_container<mydim>::bin_sort_(const UInt index, std::vector<UInt> &positions){
  // Note the scoping to avoid unnecessary storage!
  {
    std::vector<UInt> offsets{compute_offsets(index, positions)};
    for(UInt i=0; i<positions.size(); ++i){
      while(i!=offsets[i]){
        UInt next=offsets[i];
        std::swap(positions[i],positions[next]);
        std::swap(offsets[i],offsets[next]);
      }
    }
  }

  if(index>0)
    bin_sort_(index-1, positions);
}

template<UInt mydim>
std::vector<UInt> simplex_container<mydim>::compute_offsets(const UInt index, std::vector<UInt> &positions){
  std::vector<UInt> counts(num_points, 0);
  for(auto const &pos : positions)
    ++counts[simplexes[pos][index]];

  UInt offset{0};
  for (auto &count : counts){
    UInt curr{count};
    count=offset;
    offset+=curr;
  }


  std::vector<UInt> offsets;
  offsets.reserve(positions.size());
  for (auto const &pos : positions)
    offsets.push_back(counts[simplexes[pos][index]]++);
  return offsets;

}


template<UInt mydim>
void simplex_container<mydim>::check_duplicates(){
  duplicates.reserve(simplexes.size());
  // First face/edge cannot be a duplicate!
  duplicates.push_back(false);
  for(auto it=std::next(simplexes.cbegin()); it!=simplexes.cend(); ++it)
    duplicates.push_back(*std::prev(it) == *it);
}


template<UInt mydim>
void simplex_container<mydim>::store_indexes(){
  distinct_indexes.reserve(std::count(duplicates.begin(), duplicates.end(), false));
  for(UInt i=0; i<duplicates.size(); ++i)
    if(!duplicates[i])
      distinct_indexes.push_back(i);
}


template<UInt mydim>
typename simplex_container<mydim>::OutputType simplex_container<mydim>::assemble_output() const {

  std::vector<UInt> subsimplexes{assemble_subs()};
  std::vector<bool> submarkers{mark_boundary()};
  std::vector<int> neighbors{this->compute_neighbors()};

  std::vector<bool> nodesmarkers(num_points);

  for(UInt k=0; k<mydim; ++k)
    for(UInt i=0; i<submarkers.size(); ){
      nodesmarkers[subsimplexes[i+k*submarkers.size()]]=submarkers[i];
      while(i<submarkers.size() && (!submarkers[i] || nodesmarkers[subsimplexes[i+k*submarkers.size()]]))
        ++i;
    }

  return std::make_tuple(std::move(subsimplexes),std::move(submarkers),std::move(nodesmarkers),std::move(neighbors));
}

template<UInt mydim>
std::vector<UInt> simplex_container<mydim>::assemble_subs() const {
  std::vector<UInt> subsimplexes;
  subsimplexes.reserve(mydim*distinct_indexes.size());

  for(UInt j=0; j<mydim; ++j)
    for(auto const &pos : distinct_indexes)
      subsimplexes.push_back(simplexes[pos][j]);

  return subsimplexes;
}

template<UInt mydim>
std::vector<bool> simplex_container<mydim>::mark_boundary() const {
  std::vector<bool> boundarymarkers;
  boundarymarkers.reserve(distinct_indexes.size());

  std::for_each(distinct_indexes.cbegin(), std::prev(distinct_indexes.cend()), [&] (UInt i) {
    boundarymarkers.push_back(!duplicates[i+1]);
  });

  boundarymarkers.push_back(distinct_indexes.back()+1==duplicates.size() || !duplicates[distinct_indexes.back()+1]);
  return boundarymarkers;
}

template<UInt mydim>
std::vector<int> simplex_container<mydim>::compute_neighbors() const {
  std::vector<int> neighbors(simplexes.size(), -2);

  auto rep_it=duplicates.cbegin();
  simplex_t prev{simplexes.front()};
  for (auto const &curr : simplexes){
    // Note: the first simplex cannot be a duplicate!
    if (*(rep_it++)){
      neighbors[curr.i()+curr.j()*num_elements]=prev.i();
      neighbors[prev.i()+prev.j()*num_elements]=curr.i();
    }
    prev=curr;
  }
  return neighbors;
}


#endif
