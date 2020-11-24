#ifndef __AD_TREE_IMP_H__
#define __AD_TREE_IMP_H__

template<class Shape>
// Initialize header_ and data_
ADTree<Shape>::ADTree(TreeHeader<Shape> const & header): header_(header) {
  /*
   * The first element in the tree nodes vector is the head.
   * It stores the address of the tree root (i.e. the first node in the tree).
   * If it stores 0 it means that the tree is empty.
   */
  data_.reserve(header_.gettreeloc()+1);

  // Id, obj are arbitrary parameters. Remember that data_[0] is the head, not a tree node.
  Shape obj;
  Id id = std::numeric_limits<UInt>::max();
  data_.push_back(TreeNode<Shape>(id,obj));

}


//Shape is given as Element<NNODES,myDim,nDim> from mesh.h
template<class Shape>
ADTree<Shape>::ADTree(const RNumericMatrix& points, const RIntegerMatrix& triangle) {
    int ndimp = Shape::dp(); //physical dimension
    int nvertex = Shape::numVertices; //number of nodes at each Element (not total number of nodes!)

    const UInt num_nodes = points.nrows();
    const UInt num_triangle = triangle.nrows();
    // Build the tree.

    // step1: Construct TreeHeader
    std::vector<std::vector<Real> > vcoord;
    vcoord.resize(ndimp);

    for (int i = 0; i < ndimp; i++) {
      vcoord[i].resize(num_nodes);
      for(int j = 0; j < num_nodes; j++)
        vcoord[i][j] = points(j,i);
    }

    Domain<Shape> mydom(vcoord);

    //expect to have 'num_triangle' number of TreeNode
    //domain is extracted from 'points' array
    TreeHeader<Shape> myhead = createtreeheader<Shape>(num_triangle, mydom);



    //Step 2: Construct first node of adtree
    /*
     * The first element in the tree nodes vector is the head.
     * It stores the address of the tree root (i.e. the first node in the tree).
     * If it stores 0 it means that the tree is empty.
     */
    header_ = myhead;
    data_.reserve(header_.gettreeloc()+1);

    Shape obj;
    Id id = std::numeric_limits<UInt>::max();
    data_.push_back(TreeNode<Shape>(id, obj)); // Id, obj are arbitrary parameters. Remember that data_[0] is the head, not a tree node.


    // Step 3: Fill the tree: Add each element to the Treenode

    std::vector<Real> elem(nvertex*ndimp); //'elem' is a single Element, composed of vector of points
    for ( int i = 0; i < num_triangle; i++ ) {
      for (int j = 0; j < nvertex ; j++) {
        for (int l = 0; l < ndimp ; l++) {
          UInt idpt = triangle(i,j);
          elem[j*ndimp + l] =  points(idpt, l);
        }
      }
      //insert Element into tree
      this -> addtreenode(i, elem);
    } //end of loop
}


template<class Shape>
ADTree<Shape>::ADTree(SEXP Rmesh){
  int tree_loc_ = INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 3), R_DimSymbol))[0];
  int tree_lev_ = INTEGER(VECTOR_ELT(Rmesh, 11))[0];
  int ndimp_ = Shape::dp();
  int ndimt_ = 2*ndimp_;
  int nele_ = tree_loc_;
  int iava_ = tree_loc_+1;
  int iend_ = tree_loc_+1;

  std::vector<Real>  origin_;
  origin_.assign(REAL(VECTOR_ELT(Rmesh, 12)), REAL(VECTOR_ELT(Rmesh, 12))+ndimt_);
  std::vector<Real> scalingfactors_;
  scalingfactors_.assign(REAL(VECTOR_ELT(Rmesh, 13)), REAL(VECTOR_ELT(Rmesh, 13))+ndimt_);

  Domain<Shape> tree_domain(origin_, scalingfactors_);
  header_ = TreeHeader<Shape>(tree_loc_, tree_lev_, ndimp_, ndimt_, nele_, iava_, iend_, tree_domain);


  //treenode information (number of nodes = number of elements+1)
  std::vector<Id> id_;
  id_.assign(INTEGER(VECTOR_ELT(Rmesh, 14)), INTEGER(VECTOR_ELT(Rmesh, 14))+tree_loc_+1);
  std::vector<int> node_left_child_;
  node_left_child_.assign(INTEGER(VECTOR_ELT(Rmesh, 15)), INTEGER(VECTOR_ELT(Rmesh, 15))+tree_loc_+1);
  std::vector<int> node_right_child_;
  node_right_child_.assign(INTEGER(VECTOR_ELT(Rmesh, 16)), INTEGER(VECTOR_ELT(Rmesh, 16))+tree_loc_+1);
  Real* box_ = REAL(VECTOR_ELT(Rmesh, 17));

  UInt num_tree_nodes = id_.size();
  data_.reserve(num_tree_nodes);

  std::vector<Real> coord;
  coord.reserve(ndimt_);
  for (UInt i=0; i<num_tree_nodes; i++) {
    for (UInt j=0; j<ndimt_; j++) {
      coord.push_back(box_[i + num_tree_nodes*j]);
    }
    Box<Shape::dp()> box(coord);
    TreeNode<Shape> tree_node(box, id_[i], node_left_child_[i], node_right_child_[i]);
    data_.push_back(tree_node);
    coord.clear();
  }
}


template<class Shape>
int ADTree<Shape>::adtrb(Id shapeid, std::vector<Real> const & coords) {
  /*
   * We will traverse the tree in preorder fashion.
   * ipoi and ifth will store the location of the current node
   * and that of the father node.
   *
   * A value 0 for the left_link or the right_link indicates
   * a void sibling address.
   * Therefore, if left_link = right_link = 0 the node is a leaf.
   *
   * In the algorithm, when ipoi = 0 we have
   * reached a void sibling, where we can add a new node.
   */
  int nele = header_.getnele();
  int dimt = header_.getndimt();

  int iava = header_.getiava();
  int iend = header_.getiend();


  std::vector<Real> x;
  x.reserve(dimt);

  // At the start ipoi is set to the address of the root and ifth to 0.
  // If ipoi is 0, the tree is empty. Otherwise, it is the index of first true node.
  //ipoi: index to be checked next
  int ipoi = data_[0].getchild(0);
  int ifth = 0;

  /* We scale the dimension of the "bounding box" of the Shape object
   * with coordinate values given by coords.
   */

  Box<Shape::dp()> shapebox(coords);
  for(int i = 0; i < header_.getndimt(); ++i) {
    Real orig = header_.domainorig(i);
    Real scal = header_.domainscal(i);
    Real val  = (shapebox[i]-orig)*scal;
    if( (val<0.) || (val>1.) )
      // Object out of domain.
      throw(TreeDomainError<Shape>(nele+1, Shape::coordsize(), coords));
    x.push_back(val);
  }


  /* x now stores the scaled coordinates of the object's bounding box: x_i \in (0,1)
   * If Shape = Point<NDIMP> or Box<NDIMP>, the object and its bounding box coincide.
   */

  /*
   * The variable currentlev will contain the current level.
   * We start from level 0.
   */
  int currentlev = 0;
  short int edge = 0;

  while(ipoi != 0) { // finish while when ipoi == 0
    // Get the current dimension.
    int id = searchdim(currentlev, dimt);

    /*
     * We take advantage of the fact that we are only descending
     * the tree. Then we recursively multiply by 2 the coordinate.
     */

    x[id] *= 2.;
    ifth = ipoi;

    if(x[id] < 1.) {
      // Go to the left.
      edge = 0;
    } else {
      // Go to the right.
      edge = 1;
      --x[id];
    }
    // Next level.
    ++currentlev;
    ipoi = data_[ifth].getchild(edge);
    // If reached to terminal node, the getchild will always be 0 (i.e. ipoi will be 0)
    // Goal: update ifth which will be the terminal node
  }

  if( iava == iend ) {// before any deletion
    /*
     * The list of available node is empty, so we have to put the new node
     * in the yet unassigned portion of the vector storing the tree.
     */
    data_.push_back(TreeNode<Shape>(shapeid, shapebox)); // push dummy object to be changed
  }
  // else {
  //   std::vector<Real> bcoords = {shapebox[0], shapebox[1], shapebox[2], shapebox[3]};
  //   data_[iava].setcoords(bcoords);
  // }

  int neletmp = nele;
  // Add the node in the next available location.
  int ipoitmp = iava;
  // iava is the next available location.
  int iavatmp = data_[ipoitmp].getchild(0); //should return the position of left child
  ++neletmp;
  if(iavatmp == 0) {
    if( iend > header_.gettreeloc() ) {
      // Not enough space.
      throw TreeAlloc<Shape>();
    }
  }


  // Add the node in the next available location.
  ipoi = iava; // already inserted dummy object at iava
  // Set the father link to the new node.
  data_[ifth].setchild(edge, ipoi); //ifth is the terminal node index updated above


  // iava is the next available location.
  iava = data_[ipoi].getchild(0); //should return the position of left child

  ++nele;
  if(iava == 0) {
    if( iend > header_.gettreeloc() ) {
      // Not enough space.
      throw TreeAlloc<Shape>();
    }
    ++iend;
    iava = iend; //iava is updated here as next available location
  }

  /*
   * Set left_link and right link of ipoi equal to 0: this node
   * is a leaf. This is already done by the constructor of TreeNode class
   * when we put the node in the "free store".
   */

  //0 means it is empty so the children of terminal node should be 0
  data_[ipoi].setchild(0, 0);
  data_[ipoi].setchild(1, 0);
  //data_[ipoi].setfather(ifth);

  // Store back header informations.
  header_.setiend(iend);
  header_.setiava(iava);
  header_.setnele(nele);
  if(currentlev > header_.gettreelev()) {
    header_.settreelev(currentlev);
    if(currentlev > LevRuntimeError<Shape>::getmaxtreelev())
      // Current maximum number of tree levels exceeded.
      throw LevRuntimeError<Shape>();
  }

  return ipoi;
}

template<class Shape>
int ADTree<Shape>::handledomerr(Id shapeid, std::vector<Real> const & coords) {
  try {
    int iloc = adtrb(shapeid, coords);
    return iloc;
  }
  catch(TreeDomainError<Shape> de) {
    // Handle a TreeDomainError exception.
    // std::cout << "error!  " << de.getnelep1() << "-th object which is to be inserted into the tree is out of domain"
    //     << std::endl;
    // std::cout << "Coordinates" << std::endl;
    // std::cout << "-----------" << std::endl;
    // std::cout << de;
    std::exit(EXIT_FAILURE);
  }
}

template<class Shape>
int ADTree<Shape>::handletreealloc(Id shapeid, std::vector<Real> const & coords) {
  try {
    int iloc = handledomerr(shapeid, coords);
    return iloc;
  }
  catch(TreeAlloc<Shape>) {
    // Handle a TreeAlloc exception.
    // std::cout << "warning! not enough space" << std::endl;
    // std::cout << "increasing tree memory locations up to 1.5 * number of current locations..." << std::endl;
    int locv = header_.gettreeloc();
    int delta = int(locv/2);
    header_.settreeloc(locv+delta);
    if(locv == header_.gettreeloc()) {
      //std::cout << "error! no more space to add a new node" << std::endl;
      //std::exit(EXIT_FAILURE);
    }
    data_.resize(header_.gettreeloc()+1);
    int iloc = handledomerr(shapeid, coords);
    return iloc;
  }
}

template<class Shape>
int ADTree<Shape>::handleleverr(Id shapeid, std::vector<Real> const & coords) {
  try {
    int iloc = handletreealloc(shapeid, coords);
    return iloc;
  }
  catch(LevRuntimeError<Shape>) {
    // Handle a LevRuntimeError exception.
    // std::cout << "warning! maximum number of tree levels exceeded" << std::endl;
    // std::cout << "the limit is " << LevRuntimeError<Shape>::getmaxtreelev() << std::endl;
    // std::cout << "setting the new limit to" << int(LevRuntimeError<Shape>::getmaxtreelev() * 1.1) << std::endl;
    LevRuntimeError<Shape>::setmaxtreelev(int(LevRuntimeError<Shape>::getmaxtreelev() * 1.1));

    int iloc = handletreealloc(shapeid, coords);
    return iloc;
  }
}

template<class Shape>
int ADTree<Shape>::addtreenode(Id shapeid, std::vector<Real> const & coords) {
  return handleleverr(shapeid, coords);
}

template<class Shape>
void ADTree<Shape>::gettri(int const & loc, std::vector<Real> & coord, Id & id) {
  coord.clear();
  coord.reserve(Shape::dt());
  for (int i = 0; i< Shape::dt(); ++i) {
    coord.push_back(data_[loc].getcoord(i));
  }
   id = data_[loc].getid();
}

template<class Shape>
bool ADTree<Shape>::search(std::vector<Real> const & region, std::set<int> & found) const {

  
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
     tolerance = 10 * eps;

  // This function returns true if it has completed successfully, false otherwise.

  // Start preorder traversal at level 0 (root).
  int ipoi = data_[0].getchild(0);
  int ipoiNext = 0;
  int dimp = header_.getndimp();
  int dimt = header_.getndimt();

  // xl is the origin point for searching.
  std::vector<Real> xl(dimt, 0);

  std::vector<Real> box(2*dimp, 0);
  std::vector<Real> xel(2*dimp, 0);

  std::stack<std::vector<Real> > _stack;
  found.clear();

  int lev = 0;

  // Check if the region intersects the domain.
  for(int i = 0; i < dimp; ++i) {
    double orig = header_.domainorig(i);
    double scal = header_.domainscal(i);

    if(region[i] > orig+(1./scal)) { //region[i]: min of what we are searching for, orig+(1./scal): max of our domain
      return 0;
    }
    if(region[i+dimp] < orig) {  // region[i+dimp]:max of what we are searching for, orig: min of our domain
      return 0;
    }
  }

  // Rescale the region.
  // box is rescaled region
  for(int i = 0; i < dimp; ++i) {
    Real orig = header_.domainorig(i);
    Real scal = header_.domainscal(i);
    box[i] = (region[i]-orig)*scal;
    box[i+dimp] = (region[i+dimp]-orig)*scal;
  }

  /*
   * To traverse a non-empty binary tree in preorder,
   * perform the following operations recursively at each node, starting with the root node:
   * > visit the root;
   * > traverse the left subtree;
   * > traverse the right subtree.
   */
  while(ipoi != 0) {
    do {
      //xel is rescaled data_[ipoi]
      for(int i = 0; i < dimt; ++i) {
        Real orig = header_.domainorig(i);
        Real scal = header_.domainscal(i);
        xel[i] = (data_[ipoi].getcoord(i)-orig)*scal;
      }

      if(dimp == dimt) {
        //std::cout << "when dimp == dimt but in our case, there shouldn't be this case" << std::endl;
      /*
       * This function works when we have either points or boxes.
       * In the first case we have to repeat the object's coordinates.
       */
        for(int i = 0; i < dimp; ++i) xel[i+dimp] = xel[i];
      }

      // Does the element intersect box?
      int flag = 0;
      for(int i = 0; i < dimp; ++i) {
        if(box[i] > xel[i+dimp]+tolerance) { flag = 1; }
        if(box[i+dimp] < xel[i]-tolerance) { flag = 1; }
      }

      if(flag == 0) {
        found.insert(ipoi); //insert the first node that if found and then go on with the sub-tree
        /*
         * Put here all the action needed when an element is found which
         * intersects the box.
         */
      }

      // Traverse left subtree.
      ipoiNext = data_[ipoi].getchild(0);

      // Check if subtree intersects box.
      if(ipoiNext != 0) {
        int id = searchdim(lev, dimt);
        double amov = delta(lev, dimt);

        if (id < dimp) {
            if(xl[id] > box[id+dimp]) {
              ipoiNext = 0;
            }
        } else {
          if(xl[id]+amov < box[id-dimp]) {
            ipoiNext = 0;
          }
        }
      }

      /*
       * Left subtree is neither null nor 'external'.
       * Push ipoi onto the stack.
       */
      if (ipoiNext != 0) {
        std::vector<Real> stackele;
        stackele.reserve(dimt+2);
        stackele.push_back(ipoi);
        for (int i = 0; i < dimt; ++i) {
          stackele.push_back(xl[i]);
        }
        stackele.push_back(lev);
        _stack.push(stackele);
        ipoi = ipoiNext;
        ++lev;
      }

    } while(ipoiNext != 0);

    // Traverse right subtree.
    ++lev;
    ipoi = data_[ipoi].getchild(1);
    do {
      // If right_link is null we have to get the point from the stack.
      while (ipoi == 0) {
        if(_stack.empty()) { //when reached last right_link,
          return !found.empty(); //searching finished
        }
        std::vector<Real> stackele(_stack.top());
        _stack.pop();
        ipoi = data_[stackele[0]].getchild(1);
        for(int i = 0; i < dimt; ++i) {
          xl[i] = stackele[i+1];
        }
        lev = stackele[dimt+1]+1;
      }

      /*
       * Check if the subtree intersects box. Otherwise set right_link = 0,
       * pop new node from the stack and adjourn level.
       *
       * lev-1 is the level of the parent node, which directs the search.
       */
      int id = searchdim(lev-1, dimt);
      Real amov = delta(lev-1, dimt);

      // Reset xl (we have just traversed the tree to the right).
      xl[id] = xl[id]+amov;

      // Check if the subtree intersects box.
      Real x1 = xl[id];
      Real x2 = xl[id]+amov;
      if(id < dimp) {
        if(x1 > box[id+dimp]) {
          ipoi=0;
        }
      } else {
        if(x2 < box[id-dimp]) {
          ipoi=0;
        }
      }
    } while(ipoi == 0);

  } //end of while

  return !found.empty(); //if empty, return False; if not empty, return True
}

// template<class Shape>
// void ADTree<Shape>::deltreenode(int const & index) {
//   if(index < header_.getiend()) {
//     int nele = header_.getnele();
//     int ipoi = index;
//     int iava = header_.getiava();
//     int ifth = data_[ipoi].getfather();

//     --nele;
//     int lchild = data_[ipoi].getchild(0);
//     int rchild = data_[ipoi].getchild(1);
//     int whichchild;

//     if(data_[ifth].getchild(0) == ipoi) {
//       whichchild = 0;
//     } else {
//       whichchild = 1;
//     }

//     if( (lchild == 0) && (rchild == 0) ) {
//       // If the location is already a leaf, we can just release it.
//       data_[ifth].setchild(whichchild, 0);
//     } else {
//         int ipoiNext = ipoi;
//         int flag = 0;
//         int il = data_[ipoiNext].getchild(0);
//         int ir = data_[ipoiNext].getchild(1);
//         // Traverse the tree to find an empty leaf.
//         while( (il != 0) || (ir != 0) ) {
//           il = data_[ipoiNext].getchild(0);
//           ir = data_[ipoiNext].getchild(1);
//           if(il != 0) {
//             flag = 0;
//             ipoiNext = il;
//           } else if(ir != 0) {
//             flag = 1;
//             ipoiNext = ir;
//           }
//         }

//         data_[ifth].setchild(whichchild, ipoiNext);
//         int foo_f = data_[ipoiNext].getfather();
//         data_[foo_f].setchild(flag, 0);
//         data_[ipoiNext].setfather(ifth);
//         if(ipoiNext != lchild) {
//           data_[ipoiNext].setchild(0, lchild);
//           data_[lchild].setfather(ipoiNext);
//         }
//         if(ipoiNext != rchild) {
//           data_[ipoiNext].setchild(1, rchild);
//           data_[rchild].setfather(ipoiNext);
//         }
//       }

//     // Add the erased location to iava
//     data_[ipoi].setchild(0, iava);
//     data_[ipoi].setchild(1, 0);

//     // Store back header informations.
//     header_.setiava(ipoi);
//     header_.setnele(nele);
//   }
// }

template<class S>
std::ostream & operator<<(std::ostream & ostr, ADTree<S> const & myadt) {
  ostr << myadt.header_;

  return ostr;
}


#endif /* ADTREE_IMP_H_ */
