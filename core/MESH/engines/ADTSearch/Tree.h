#ifndef __TREE_H__
#define __TREE_H__

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <queue>
#include <functional>

#include <iostream>


// An implementation of a binary tree structure

// forward declaration
template <typename T> class Node;
// define a node_ptr as a shared pointer to a node
template <typename T> using node_ptr = std::shared_ptr<Node<T>>;

enum LinkDirection {LEFT, RIGHT};

// T is the type of the object stored in the node, N the number of children a node can have.
// Use N = 2 for binary trees
template <typename T> class Node {
private:
  
  T data_;                               // the actual stored data
  std::array<node_ptr<T>, 2> children_;  // children_[0] is the left child. children_[1] is the right child
  unsigned int key_;                     // unique node key

public:
  // constructor
  Node(T data, unsigned int key) : data_(data), key_(key) {} ;

  // add child at first available position. Returns nullptr if there is no room
  node_ptr<T> addChild(const T& data, unsigned int key);
  // add child at node at the specified position. Returns nullptr if the index is already in use
  node_ptr<T> addChild(const T& data, unsigned int key, LinkDirection index);

  // a method to check if this node has no children (leaf node)
  bool isLeaf() const;

  // getters
  T getData() const { return data_; }
  std::array<node_ptr<T>, 2> getChildren() const { return children_; }
  unsigned int getKey() const { return key_; }
};

template <typename T> node_ptr<T> Node<T>::addChild(const T& data, unsigned int key) {
  for (size_t i = 0; i < children_.size(); ++i){
    if(children_[i] == nullptr){
      // add node at first available position
      children_[i] = std::make_shared<Node<T>>(data, key);
      
      return children_[i];
    }
  }
  return nullptr;
}

template <typename T> node_ptr<T> Node<T>::addChild(const T& data, unsigned int key, LinkDirection index) {
  // link already in use
  if(children_[index] != nullptr)
    return nullptr;
  
  // add node as child
  children_[index] = std::make_shared<Node<T>>(data, key);
  return children_[index];
}

template <typename T> bool Node<T>::isLeaf() const {
  // check if all pointers are null. If not, this is not a leaf
  for(node_ptr<T> child : children_){
    if(child != nullptr)
      return false;
  }
  return true;
}

// a binary tree storing objects of type T.
template <typename T>
class Tree {

private:

  // pointer to root (from the root you have access to any other node)
  node_ptr<T> root_;
  // counter to keep track of how many nodes are stored in the tree. It also provides
  // a way for unique indexing nodes by ID
  unsigned int numNodes = 0;

  // auxiliary data structure mapping node indexes to pointers to nodes
  std::unordered_map<unsigned int, node_ptr<T> > nodeTable;
  
public:
  // an empty tree
  Tree() = default;
  // initialize the tree with its root node
  Tree(const T& root) {
    node_ptr<T> rootNode = std::make_shared<Node<T>>(root, 0);
    root_ = rootNode;

    // store reference to root in nodeTable
    nodeTable[numNodes] = rootNode;
    numNodes++;
  }

  // insert node at first available position
  void insert(const T& data);
  // insert node as child given the ID of the father. Do nothing if there is already a child in that direction.
  // This routine is usefull if the tree has to be built progressively using some external criterion
  bool insert(const T& data, unsigned int ID, LinkDirection direction);
  
  // get node given its ID
  node_ptr<T> getNode(unsigned int ID) const { return nodeTable.at(ID); } ;

  unsigned int getNumberOfNodes() const { return numNodes; }
};

// insert a node in the first available position
template <typename T>
void Tree<T>::insert(const T& data) {
  // perform a level-order traversal to find first available position
  std::queue<unsigned int> queue;

  // insert root key
  queue.push(0);

  while(!queue.empty()){
    // get pointer to node
    unsigned int key = queue.front();
    queue.pop();
    node_ptr<T> node = nodeTable[key];
    
    // this evaluates true if the insertion happened
    node_ptr<T> newNode = node->addChild(data, numNodes);
    if(newNode != nullptr){
      nodeTable[numNodes] = newNode;
      numNodes++;
      return;
    }else{
      // in case the node is already full add its children ID to the queue for later processing
      for(node_ptr<T> child : node->getChildren()){
	queue.push(child->getKey());
      }
    }
  }
  return;
}

template <typename T>
bool Tree<T>::insert(const T &data, unsigned int ID, LinkDirection direction) {
  // take father node using nodeTable
  node_ptr<T> father = nodeTable.at(ID);

  // add child to father at given direction
  node_ptr<T> child = father->addChild(data, numNodes, direction);
  if (child != nullptr){ // insertion is ok
    // add child to nodeTable
    nodeTable[numNodes] = child;
    numNodes++;
    return true;
  }  
  return false;
}

#endif // __TREE_H__
