#ifndef __TREE_HEADER_IMP_H__
#define __TREE_HEADER_IMP_H__


template<class T>
TreeHeader<T>::TreeHeader(int const & ntree, Domain<T> const & d):
	tree_loc_(ntree), tree_lev_(0), ndimp_(T::dp()), ndimt_(T::dt()), nele_(0), iava_(1), iend_(1), tree_domain_(d) {

	std::vector<TreeNode<T> > foo;
	if(foo.max_size() < unsigned(tree_loc_ + 1))
		/* If there is no enough space to store the requested nodes and
		 * the tree head, a LocLengthError exception is thrown.
		 */
		throw(LocLengthError<T>(foo.max_size(), tree_loc_));
}

template<class T>
void TreeHeader<T>::stml(int const & nt) {
	std::vector<TreeNode<T> > foo;
	if(foo.max_size() < unsigned(nt + 1))
		/* If there is no enough space to store the requested nodes and
		 * the tree head, a LocLengthError exception is thrown.
		 */
		throw(LocLengthError<T>(foo.max_size(), nt));
	tree_loc_ = nt;
}

template<class T>
std::ostream & operator<<(std::ostream & ostr, TreeHeader<T> const & head) {
	ostr << "General informations about the tree" << std::endl;
	ostr << "----------------------------------" << std::endl;
	ostr << "Tree memory locations: " << head.tree_loc_ << std::endl;
	ostr << "Number of tree levels: " << head.tree_lev_ << std::endl;
	ostr << "Number of physical space dimension: " << head.ndimp_ << std::endl;
	ostr << "Number of dimensions used for the search: " << head.ndimt_ << std::endl;
	ostr << "Number of logical locations currently used in the tree: " << head.nele_ << std::endl;
	ostr << head.tree_domain_;
	ostr << "----------------------------------" << std::endl;

	return ostr;
}

template<class T>
TreeHeader<T> createtreeheader(int const & nt, Domain<T> const & d) {
	try {
		TreeHeader<T> hd(nt, d);
		return hd;
	}
	catch(LocLengthError<T> lo) {
		// std::cout << std::endl << std::endl;
		// std::cout << "warning!	createtreeheader : max dimension exceeded" << std::endl;
		// std::cout << "the limit is " << lo.getmaxtreeloc()-1
		// 		<< " while needed at least " << lo.gettreeloc() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

template<class T>
void TreeHeader<T>::settreeloc(int const & nt) {
	try {
		stml(nt);
	}
	catch(LocLengthError<T> lo) {
		// std::cout << std::endl << std::endl;
		// std::cout << "warning!	settreeloc : max dimension exceeded" << std::endl;
		// std::cout << "the limit is " << lo.getmaxtreeloc()-1
		// 		<< " while requested " << lo.gettreeloc() << std::endl;
		// std::cout << "increasing tree memory locations up to the limit" << std::endl;
		stml(lo.getmaxtreeloc()-1);
	}
}

#endif /* TREE_HEADER_IMP_H_ */
