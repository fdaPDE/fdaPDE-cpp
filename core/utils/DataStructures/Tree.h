#ifndef __TREE_H__
#define __TREE_H__

#include <memory>
#include <unordered_map>
#include <queue>

namespace fdaPDE{
namespace core{
  
  // An implementation of a binary tree data structure

  // forward declaration
  template <typename T> class Node;
  // define a node_ptr as a shared pointer to a node
  template <typename T> using node_ptr = std::shared_ptr<Node<T>>;

  enum LinkDirection {LEFT, RIGHT};

  // a tree node holding a general type T object
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
    T getData()                              const { return data_;     }
    std::array<node_ptr<T>, 2> getChildren() const { return children_; }
    unsigned int getKey()                    const { return key_;      }
  };

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

    // insert node at first available position, returns the pointer to the inserted node
    node_ptr<T> insert(const T& data);
    // insert node as child given the ID of the father. Do nothing and return nullptr if there is already a child in that direction.
    // This routine is usefull if the tree has to be built progressively using some external criterion
    node_ptr<T> insert(const T& data, unsigned int ID, LinkDirection direction);
  
    // getters
    node_ptr<T> getNode(unsigned int ID) const { return nodeTable.at(ID); }
    unsigned int getNumberOfNodes()      const { return numNodes;         }
  };

#include "Tree.tpp"
}}
  
#endif // __TREE_H__
