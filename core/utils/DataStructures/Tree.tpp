// Node

// add node containing data as child of "key" at first available position. Do nothing if no child can be added
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

// add node containing data as child of "key" at the given position (LEFT or RIGHT child). Do nothing if link
// already busy
template <typename T> node_ptr<T> Node<T>::addChild(const T& data, unsigned int key, LinkDirection index) {
  // link already in use
  if(children_[index] != nullptr)
    return nullptr;
  
  // add node as child
  children_[index] = std::make_shared<Node<T>>(data, key);
  return children_[index];
}

// check if this node is a leaf
template <typename T> bool Node<T>::isLeaf() const {
  // check if all pointers are null. If not, this is not a leaf
  for(node_ptr<T> child : children_){
    if(child != nullptr)
      return false;
  }
  return true;
}

// Tree

// insert a node at first available position
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

// add node as child of ID following the specified direction (LEFT or RIGHT child). Returns false if
// the operation cannot be performed.
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
