#ifndef __TREE_NODE_IMP_H__
#define __TREE_NODE_IMP_H__


template<class T>
void TreeNode<T>::print(std::ostream & out) const
{
	//Be careful! this id is not the id of the Treenode but id of Triangle
	out << "------------------------------" << std::endl;
	out << "Shape Id: --" << id_ << "--" <<std::endl;
	//out << "Father:  " << father_ << std::endl;
	out << "Left Children:  " << children_[0] << std::endl;
	out << "Right Children:  " << children_[1] << std::endl;
	out << "Box: ";
	box_.print(out);
	out<<std::endl;
}

#endif /* TREENODE_IMP_H_ */
