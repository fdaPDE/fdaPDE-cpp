#ifndef __TREE_SEARCH_H__
#define __TREE_SEARCH_H__

#include <memory>
#include <unordered_map>
#include <queue>

// T is the type of the object stored in the node, N the number of children a node can have.
// Use N = 2 for binary trees
template <typename T, unsigned int N> class Node {
private:
  // the actual stored data
  T data_;
  // for binary trees, childs_[0] is always the left child. childs_[1] the right one
  std::array<std::shared_ptr<Node>, N> children_;
  // unique node key
  unsigned int key_;
  
public:
  Node(T data, unsigned int key) : data_(data), key_(key) {} ;

  // add child at node at the specified position. Returns false if the index is out of range
  bool addChild(const Node& child, unsigned int index);

  // a method to check if this node has no children (leaf node)
  bool isLeaf() const;
  
  T getData() const { return data_; }
  std::array<std::shared_ptr<Node>, N> getChildren() const { return children_; }
  
};



// be careful, if the index is already in use the node will be substituted
template <typename T, unsigned int N> bool Node<T, N>::addChild(const Node<T,N>& child, unsigned int index) {
  // index out of range
  if(index > N)
    return false;
  
  // add the node in the list of children
  children_[index] = std::make_shared<Node<T, N>>(child);
  return true;
}

template <typename T, unsigned int N> bool Node<T, N>::isLeaf() const {

  // check if all pointers are null. If not this is not a leaf
  for(std::shared_ptr<Node<T,N>> child : children_){
    if(!child != nullptr)
      return false;
  }

  return true;
}

// a tree storing objects of type T. N indicates the degree of the tree
// Use N = 2 for binary trees
template <typename T, unsigned int N>
class Tree {

private:

  // pointer to root (from the root you have access to any other node)
  std::shared_ptr<Node<T,N>> root_;
  // counter to keep track of how many nodes are stored in the tree. It also provides
  // a way for unique indexing nodes by ID
  unsigned int numNodes = 0;

  // auxiliary data structure mapping node indexes to pointers to nodes
  std::unordered_map<unsigned int, std::shared_ptr<Node<T,N>> > nodeTable;
  
public:
  // an empty tree
  Tree() = default;
  // initialize the tree with its root node
  Tree(const T& root) {
    std::shared_ptr<Node<T,N>> rootNode = std::make_shared<Node<T,N>>(root, 0);
    root_ = rootNode;

    // store reference to root in nodeTable
    nodeTable[numNodes] = rootNode;
    numNodes++;
  }

  void insert(const T& node);

  // insert node at first free position.
  // insert node at location if possible, first free child
  // insert node at location if possible, child position choosen
  // use some custom policy to navigate the tree based on T and insert at first free position
};

// insert a node in the first available position
template <typename T, unsigned int N>
void Tree<T,N>::insert(const T&node) {
  // perform a level-order traversal to find first available position
  std::queue<unsigned int> queue;

  // insert root key
  queue.push(0);

  while(!queue.empty()){
    // get pointer to node
    unsigned int key = queue.front();
    queue.pop();
    std::shared_ptr<Node<T,N>> node = nodeTable[key];

    // create a method in node to add a child at first available position.
    // call here add child, if it returns false, add childs of node to the queue
    // executing this recursively performs a level-order traversal of the tree.
  }
}


#endif
