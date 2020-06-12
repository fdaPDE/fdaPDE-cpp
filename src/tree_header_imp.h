#ifndef TREE_HEADER_IMP_H_
#define TREE_HEADER_IMP_H_


template<UInt ndim>
TreeHeader<ndim>::TreeHeader(int const & ntree, Domain<ndim> const & d) :
	tree_loc_(ntree), tree_domain_(d) {

	std::vector<TreeNode<ndim> > foo;
	if(foo.max_size() < unsigned(tree_loc_ + 1))
		/* If there is not enough space to store the requested nodes and
		 * the tree head, a LocLengthError exception is thrown.
		 */
		throw(LocLengthError<T>(foo.max_size(), tree_loc_));
}

template<UInt ndim>
void TreeHeader<ndim>::stml(int const & nt) {
	std::vector<TreeNode<ndim> > foo;
	if(foo.max_size() < unsigned(nt + 1))
		/* If there is no enough space to store the requested nodes and
		 * the tree head, a LocLengthError exception is thrown.
		 */
		throw(LocLengthError<T>(foo.max_size(), nt));
	tree_loc_ = nt;
}

template<UInt NDIM>
std::ostream & operator<<(std::ostream & ostr, TreeHeader<NDIM> const & head) {
	ostr << "General informations about the tree" << std::endl;
	ostr << "----------------------------------" << std::endl;
	ostr << "Tree memory locations: " << head.tree_loc_ << std::endl;
	ostr << "Number of tree levels: " << head.tree_lev_ << std::endl;
	ostr << "Number of physical space dimension: " << NDIM << std::endl;
	ostr << "Number of dimensions used for the search: " << 2*NDIM << std::endl;
	ostr << "Number of logical locations currently used in the tree: " << head.nele_ << std::endl;
	ostr << head.tree_domain_;
	ostr << "----------------------------------" << std::endl;

	return ostr;
}

template<UInt NDIM>
TreeHeader<NDIM> createtreeheader(int const & nt, Domain<NDIM> const & d) {
	try {
		TreeHeader<NDIM> hd(nt, d);
		return hd;
	}
	catch(LocLengthError<T> lo) {
		std::cout << std::endl << std::endl;
		std::cout << "warning!	createtreeheader : max dimension exceeded" << std::endl;
		std::cout << "the limit is " << lo.getmaxtreeloc()-1
				<< " while needed at least " << lo.gettreeloc() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

template<UInt ndim>
void TreeHeader<ndim>::settreeloc(int const & nt) {
	try {
		stml(nt);
	}
	catch(LocLengthError<T> lo) {
		std::cout << std::endl << std::endl;
		std::cout << "warning!	settreeloc : max dimension exceeded" << std::endl;
		std::cout << "the limit is " << lo.getmaxtreeloc()-1
				<< " while requested " << lo.gettreeloc() << std::endl;
		std::cout << "increasing tree memory locations up to the limit" << std::endl;
		stml(lo.getmaxtreeloc()-1);
	}
}

#endif /* TREE_HEADER_IMP_H_ */
