#ifndef __MBTree_H__
#define __MBTree_H__

#include "Tree.h"
#include <memory>
#include <vector>

namespace fdaPDE{
namespace core{

  // A multi tree is a case of DAG which encodes the idea of "multiple overlapping
  // trees" A multi binary tree (MBT) is the case in which each tree in the multi
  // tree is a binary tree. This information allow us to take different binary
  // trees which have at least one node, as well as the whole subtree rooted in
  // such node, in common and compress their storage in a unique structure. Given
  // that each subtree in the MBT is binary, from an MBT we can extract each of
  // the original trees.

  /* For example the following 3 binary trees
  
               N02            |             N12            |             N22          
        N01           N11     |      N11           N21     |      N21           N31   
     N00   N10     N10   N20  |   N10   N20     N20   N30  |   N20   N30     N30   N40
                              |                            |
                A             |              B             |              C

     are compressed in the following single MBT

               N02            N12           N22
        N01           N11            N21           N31   
     N00   N10     N10   N20      N20   N30     N30   N40
     |---------------------|
                A
		   |----------------------|
		               B
			          |---------------------|
				             C
  */
  // used in the implementation of Splines (for example the above represents an MBT for a set of 3 cubic B-splines)

  template <typename T>
  class MBTree {
  private:
    // vector of roots (store ID of roots)
    std::vector<unsigned int> roots_;
    
    // auxiliary data structure mapping node indexes to pointers to nodes
    std::unordered_map<unsigned int, node_ptr<T>> nodeTable_;

  public:
    // constructor initializes an empty MBTree
    MBTree() = default;
    
    // insert root in the tree
    void insertRoot(const T& data, unsigned int rootID);

    // insert data at first available location under root tree with identifier ID
    node_ptr<T> insert(unsigned int root, const T& data, unsigned int ID);
    // insert data as child node of fatherID with identifier nodeID along the specified link direction. If the node already exists
    // no node is created and just the proper pointers are updated.
    node_ptr<T> insert(unsigned int fatherID, const T& data, unsigned int nodeID, LinkDirection direction);
    
    // extract the binary tree rooted under root as a View
    const Tree<T> treeView(unsigned int rootID) const;
  };

  template <typename T>
  void MBTree<T>::insertRoot(const T& data, unsigned int rootID){
    if(nodeTable_.count(rootID)) return; // node with ID rootID already present, do nothing

    // create node for root
    node_ptr<T> rootNode = std::make_shared<Node<T>>(data, rootID);
    roots_.push_back(rootID);         // store root ID in roots_vector
    
    // store reference to root in nodeTable_
    nodeTable_[rootID] = rootNode;
    return;
  }

  // insert data at first available location under root tree with identifier ID
  template <typename T>
  node_ptr<T> MBTree<T>::insert(unsigned int root, const T& data, unsigned int ID){
    // perform a level-order traversal to find first available position
    std::queue<unsigned int> queue;
    queue.push(root);    // insert root key

    while(!queue.empty()){
      // get pointer to node
      unsigned int key = queue.front();
      queue.pop();
      node_ptr<T> node = nodeTable_[key];
    
      // this evaluates true if the insertion happened
      node_ptr<T> newNode = node->addChild(data, ID);
      if(newNode != nullptr){
	nodeTable_[ID] = newNode;
	return newNode;
      }else{
	// in case the node is already full add its children ID to the queue for later processing
	for(node_ptr<T> child : node->getChildren()){
	  queue.push(child->getKey());
	}
      }
    }
    return nullptr;
  }

  // insert data as child given the ID of the father and a link direction. Do nothing if there is already a node in that direction.
  template <typename T>
  node_ptr<T> MBTree<T>::insert(unsigned int fatherID, const T& data, unsigned int nodeID, LinkDirection direction) {
    // take father node using nodeTable
    node_ptr<T> father = nodeTable_.at(fatherID);

    if(nodeTable_.count(nodeID)){ // node already present, update pointers without actual creating any node. data argument is discarded
      father->set_ptr(nodeTable_[nodeID], direction);
      return nodeTable_[nodeID];
    }else{ // create node and update MBTree
      node_ptr<T> child = father->addChild(data, nodeID, direction); // add child to father at given direction
      if (child != nullptr){ // insertion is ok
	// add child to nodeTable_
	nodeTable_[nodeID] = child;
	return child;
      }  
      return nullptr;
    }
  }

  // extract the binary tree rooted under root. This returns a view over the MBTree.
  template <typename T>
  const Tree<T> MBTree<T>::treeView(unsigned int rootID) const{
    return Tree<T>(nodeTable_.at(rootID));
  }
  
}}

#endif // __MBTree_H__
