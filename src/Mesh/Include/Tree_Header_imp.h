#ifndef __TREE_HEADER_IMP_H__
#define __TREE_HEADER_IMP_H__


template<class T>
TreeHeader<T>::TreeHeader(int const & ntree, Domain<T> const & d):
	tree_loc_(ntree), tree_lev_(0), ndimp_(T::dp()), ndimt_(T::dt()), nele_(0), iava_(1), iend_(1), tree_domain_(d) {}

template<class T>
void TreeHeader<T>::stml(int const & nt) {
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
		TreeHeader<T> hd(nt, d);
		return hd;
}

template<class T>
void TreeHeader<T>::settreeloc(int const & nt) {
		stml(nt);
}

#endif /* TREE_HEADER_IMP_H_ */
