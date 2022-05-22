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
node_ptr<T> Tree<T>::insert(const T& data) {
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

// add node as child of ID following the specified direction (LEFT or RIGHT child). Returns false if
// the operation cannot be performed.
template <typename T>
node_ptr<T> Tree<T>::insert(const T &data, unsigned int ID, LinkDirection direction) {
  // take father node using nodeTable
  node_ptr<T> father = nodeTable.at(ID);

  // add child to father at given direction
  node_ptr<T> child = father->addChild(data, numNodes, direction);
  if (child != nullptr){ // insertion is ok
    // add child to nodeTable
    nodeTable[numNodes] = child;
    numNodes++;
    return child;
  }  
  return nullptr;
}

template <typename T>
template <typename F>
void Tree<T>::DFS(const F& f) const {
  std::stack<node_ptr<T>> stack{};        // use a stack to assist the search
  std::set<unsigned int> visitedNodes{};  // set containing the IDs of visited nodes
    
  // push reference to root on stack
  stack.push(root_);

  // start visit
  while(!stack.empty()){ // repeat until stack is not empty
    node_ptr<T> currentNode = stack.top();
    stack.pop();
    if(visitedNodes.find(currentNode->getKey()) == visitedNodes.end()){ // node still not marked as visited
      visitedNodes.insert(currentNode->getKey()); // mark this node as visited

      // execute functor on data contained in this node
      f(currentNode->getData());

      // add child nodes to stack for later processing
      for(const auto& childNode : currentNode->getChildren()){
	if(childNode != nullptr)
	  stack.push(childNode);
      }
    }
  }
  return;
}
