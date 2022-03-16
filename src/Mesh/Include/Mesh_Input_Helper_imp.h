#ifndef __MESH_INPUT_HELPER_IMP_H__
#define __MESH_INPUT_HELPER_IMP_H__


template<UInt mydim>
template<std::size_t SIZE>
void simplex_container<mydim>::fill_container(const std::array<UInt, SIZE>& ORDERING){
 static_assert(SIZE==mydim*(mydim+1) || (mydim==2 && SIZE==12),
        "ERROR! ORDERING SIZE SHOULD BE EQUAL TO 2X THE NUMBER OF EDGES OR 3X THE NUMBER OF FACES! See: mesh_input_helper_imp.h");

 const UInt num_elements=elements.nrows();

 simplexes.reserve(num_elements*ORDERING.size()/mydim);

 {
   std::array<UInt,mydim> curr;
   for(UInt i=0; i<num_elements; ++i){
     for(UInt j=0; j<ORDERING.size()/mydim; ++j){
        for(UInt k=0; k<mydim; ++k)
          curr[k]=elements(i, ORDERING[mydim*j+k]);
        std::sort(curr.begin(), curr.end());
        simplexes.emplace_back(i,j,curr);
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
  
  const UInt num_points=nodes.nrows();

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
void simplex_container<mydim>::assemble_subs(SEXP Routput, UInt index) const {
  
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, distinct_indexes.size(), mydim));
  RIntegerMatrix subsimplexes(VECTOR_ELT(Routput, index));

  for(UInt j=0; j<mydim; ++j)
    for(UInt i=0; i<distinct_indexes.size(); ++i)
      subsimplexes(i,j) = simplexes[distinct_indexes[i]][j] + 1;

}

template<UInt mydim>
void simplex_container<mydim>::mark_boundary(SEXP Routput, UInt index) const {
  
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(LGLSXP, distinct_indexes.size(), 1));
  RIntegerMatrix boundarymarkers(VECTOR_ELT(Routput, index));

  for(UInt i=0; i<distinct_indexes.size()-1; ++i)
    boundarymarkers[i] = !duplicates[distinct_indexes[i]+1];

  //Special attention for the last simplex!
  boundarymarkers[distinct_indexes.size()-1] = distinct_indexes.back()+1==duplicates.size() || !duplicates[distinct_indexes.back()+1];
}


template<UInt mydim>
void simplex_container<mydim>::compute_neighbors(SEXP Routput, UInt index) const {
  
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, simplexes.size()/(mydim+1), mydim+1));
  RIntegerMatrix neighbors(VECTOR_ELT(Routput, index));

  for (UInt i=0; i<simplexes.size(); ++i)
    neighbors[i]=-1;

  auto rep_it=duplicates.cbegin();
  simplex_t prev{simplexes.front()};
  for (auto const &curr : simplexes){
    // Note: the first simplex cannot be a duplicate!
    if (*(rep_it++)){
      neighbors(curr.i(), curr.j()) = prev.i()+1;
      neighbors(prev.i(), prev.j()) = curr.i()+1;
    }
    prev=curr;
  }
}

template<UInt mydim>
void simplex_container<mydim>::order2extend(SEXP Routput, UInt index) const {
  static_assert(mydim==2, 
    "ERROR! ORDER 2 EXTENSIONS IS INTENDED FOR EDGE CONTAINERS ONLY! See mesh_input_helper_imp");
  
  const UInt num_extra_nodes = (isTriangleContainer) ? 3 : 6;

  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, simplexes.size()/num_extra_nodes, num_extra_nodes));
  RIntegerMatrix edges_extended(VECTOR_ELT(Routput, index));

  {
    UInt offset{nodes.nrows()};
    UInt pos=0;
    for(auto const &curr : simplexes){
      offset += !duplicates[pos++];
      edges_extended(curr.i(), curr.j()) = offset;
    }
  }
}


template<UInt mydim>
std::vector<UInt> simplex_container<mydim>::how_many_neighbors(UInt index) const{
    static_assert(mydim==1,
                  "ERROR! how_many_neighbors IS INTENDED FOR POINTS CONTAINERS ONLY! See mesh_input_helper_imp");

    std::vector<UInt> res;

     UInt curr_idx = index;
     while( curr_idx < simplexes.size() && simplexes[index] == simplexes[curr_idx] )
     {
        //res is never empty!
        res.push_back( curr_idx );
        ++curr_idx;
     }
    return res;
}

//Specialization is needed.
template<>
void simplex_container<1>::compute_neighbors(SEXP Routput, UInt index) const {

    //Matrix storing how many neighbors each element has, according to the side
    SET_VECTOR_ELT(Routput,index,Rf_allocMatrix(INTSXP,elements.nrows(),2));
    RIntegerMatrix lengths(VECTOR_ELT(Routput,index));

    for(UInt i : this->distinct_indexes){
        //vector of all indexes of simplexes sharing  the same node (i.e. simplexes[distinc_indexes[i]].node[0])
        //nb) currents is never empty.
        std::vector<UInt> currents = this->how_many_neighbors(i);
        UInt tot_neighbors =  currents.size() - 1;
        for( UInt idx : currents ){
            lengths( simplexes[idx].i() , simplexes[idx].j() ) = tot_neighbors;
        }
    }


    //
    //     * - - e1 - - * - - e2 - - *        edges
    //     0           1  0          1 nodes (local numbering)
    //    (1)        (0)  (1)       (0) sides (local numbering)
    //each row of neighbors matrix contains two lists of edges
    //j = 0 there is the "left" list
    //j = 1 there is the "right" list

    SET_VECTOR_ELT(Routput,index+1,Rf_allocMatrix(VECSXP,elements.nrows(),2));
    for(UInt i=0; i<elements.nrows() * 2; ++i){
        SET_VECTOR_ELT( VECTOR_ELT(Routput,index+1) ,i,Rf_allocMatrix(INTSXP,1,lengths[i]));
    }

    RIntMatrixMatrix neighbors(VECTOR_ELT(Routput,index+1));

    //Filling neighbors matrix
    for(UInt i : this->distinct_indexes) {
        std::vector <UInt> currents = this->how_many_neighbors(i);

        for (UInt j: currents) {
            UInt pos = 0;
            for (UInt k: currents) {
                if (k != j) {
                    //Indexes in R starts from 1, in C++ from 0, needed transformations!
                    neighbors(simplexes[j].i(), simplexes[j].j())[pos] = simplexes[k].i() + 1;
                    ++pos;
                }
            }
        }
    }

}
#endif
