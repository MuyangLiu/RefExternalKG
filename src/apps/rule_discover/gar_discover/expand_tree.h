#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPAND_TREE_H
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPAND_TREE_H

#include <string>
#include <sstream>
#include <vector>

namespace grape {

namespace _gar_discover {

template<typename GraphPatternType,
         typename    DataGraphType>
class ExpandTreeNode;

template<typename GraphPatternType,
         typename    DataGraphType>
class ExpandTreeLevel;

template <typename GraphPatternType,
         typename    DataGraphType>
std::string& operator<<(
    std::string& out_string,
    const ExpandTreeNode<GraphPatternType,
                            DataGraphType>& expand_tree_node) {
  out_string = std::move(out_string) + " <expand_tree_node";

  out_string = std::move(out_string)  
             + " " + GUNDAM::ToString(expand_tree_node.id());

  out_string = std::move(out_string) + " ";

  std::string str;

  str << expand_tree_node.const_pattern();

  out_string = std::move(out_string) + str;

  // serilize new vertexes
  out_string = std::move(out_string) + " <new_vertexes";
  for (auto vertex_id_cit  = expand_tree_node.NewVertexesIdCBegin();
            vertex_id_cit != expand_tree_node.NewVertexesIdCEnd();
            vertex_id_cit++){
    out_string = std::move(out_string) + " " + std::to_string(*vertex_id_cit);
  }
  out_string = std::move(out_string) + " >"; // new_vertexes

  // serilize possible rhs literal set
  out_string = std::move(out_string) + " <rhs_literals";
  for (auto literal_cit  = expand_tree_node.RhsLiteralsBegin();
            literal_cit != expand_tree_node.RhsLiteralsEnd();
            literal_cit++){
    out_string = std::move(out_string) + " ";
    out_string << *literal_cit;
  }
  out_string = std::move(out_string) + " >"; // rhs_literals

  // serilize all lhs literal set
  out_string = std::move(out_string) + " <lhs_literals";
  for (auto literal_cit  = expand_tree_node.LhsLiteralsBegin();
            literal_cit != expand_tree_node.LhsLiteralsEnd();
            literal_cit++){
    out_string = std::move(out_string) + " ";
    out_string << *literal_cit;
  }
  out_string = std::move(out_string) + " >"; // lhs_literals

  out_string = std::move(out_string) + " >"; // expand_tree_node

  return out_string;
}

template <typename GraphPatternType,
          typename    DataGraphType>
std::string& operator>>(
    std::string& in_string,
    ExpandTreeNode<GraphPatternType, 
                      DataGraphType>& expand_tree_node) {
  using namespace GUNDAM;

  using DataGraphLiteralType = gar::LiteralInfo<GraphPatternType,
                                                   DataGraphType>;

  using PatternVertexIDType = typename GraphPatternType
                                           ::VertexType
                                               ::IDType;

  std::stringstream ss;
  ss << in_string;

  std::string str;

  ss >> str;
  assert(str == "<expand_tree_node");

  // deserilize expand node id
  ss >> str;
  int expand_node_id = GUNDAM::StringToDataType<int>(str);

  // deserilize graph pattern
  GraphPatternType graph_pattern;
  getline(ss, str);
  str >> graph_pattern;
  ss.clear();
  ss << str;

  // deserilize new vertexes id set
  std::vector<PatternVertexIDType> new_vertexes_id;
  ss >> str;
  assert(str == "<new_vertexes");
  while (ss >> str) {
    if (str == ">") {
      break;
    }
    new_vertexes_id.emplace_back(
        GUNDAM::StringToDataType<PatternVertexIDType>(str));
  }
  assert(new_vertexes_id.size() == 0
      || new_vertexes_id.size() == 1);

  // deserilize possible rhs literal set
  std::vector<DataGraphLiteralType> rhs_literals;
  ss >> str;
  assert(str == "<rhs_literals");

  while (ss.peek() == ' ') {
    char c = ss.get();
  }

  while (ss.peek() != '>') {
    DataGraphLiteralType literal;
    getline(ss, str);
    str >> literal;
    ss.clear();
    ss << str;
    rhs_literals.emplace_back(literal);
    while (ss.peek() == ' ') {
      char c = ss.get();
    }
  }

  ss >> str;
  assert(str == ">");  // rhs_literals

  // deserilize all lhs literal set
  std::vector<DataGraphLiteralType> lhs_literals;
  ss >> str;
  assert(str == "<lhs_literals");

  while (ss.peek() == ' ') {
    char c = ss.get();
  }

  while (ss.peek() != '>') {
    DataGraphLiteralType literal;
    getline(ss, str);
    str >> literal;
    ss.clear();
    ss << str;
    lhs_literals.emplace_back(literal);
    while (ss.peek() == ' ') {
      char c = ss.get();
    }
  }

  ss >> str;
  assert(str == ">");  // lhs_literals

  ss >> str;
  assert(str == ">");  // expand_tree_node

  getline(ss, in_string);
  if (ss.fail()) 
    in_string.clear();

  if (new_vertexes_id.empty()){
    // does not have new vertexes
    expand_tree_node = ExpandTreeNode<GraphPatternType, 
                                         DataGraphType>(expand_node_id,
                                                        graph_pattern,
                                                        rhs_literals,
                                                        lhs_literals);
    return in_string;
  }
  assert(new_vertexes_id.size() == 1);
  // has new vertexes
  expand_tree_node = ExpandTreeNode<GraphPatternType, 
                                       DataGraphType>(expand_node_id,
                                                      graph_pattern,
                                                     *new_vertexes_id.begin(),
                                                      rhs_literals,
                                                      lhs_literals);
  return in_string;
}

// the node for the expand tree
// the expand tree is different from the generation tree, which does not
// hold the literal tree in each node
// 
// the most essential data hold in it is the graph pattern, the other
// data are just used for pruning and acceleration
template<typename _GraphPatternType,
         typename     DataGraphType>
class ExpandTreeNode{
 public:
  using GraphPatternType = _GraphPatternType;

 private:
  friend ExpandTreeLevel<_GraphPatternType,
                             DataGraphType>;

  using DataGraphLiteralType = gar::LiteralInfo<_GraphPatternType,
                                                    DataGraphType>;

  using PatternVertexIDType = typename GraphPatternType
                                           ::VertexType
                                               ::IDType;

  inline void set_id(int id){
    this->id_ = id;
    return;
  }

 public:
  ExpandTreeNode() = default;

  // construct method with no new vertex
  ExpandTreeNode(int id,
                 const GraphPatternType& graph_pattern,
                 const std::vector<DataGraphLiteralType>& rhs_literals,
                 const std::vector<DataGraphLiteralType>& lhs_literals)
                :id_(id),
                 graph_pattern_(graph_pattern),
                 new_vertexes_id_(),
                 rhs_literals_(rhs_literals),
                 lhs_literals_(lhs_literals){
    return;
  }

  // construct method with one new vertex
  ExpandTreeNode(int id,
                 const GraphPatternType& graph_pattern,
                 const PatternVertexIDType& new_vertex_id,
                 const std::vector<DataGraphLiteralType>& rhs_literals,
                 const std::vector<DataGraphLiteralType>& lhs_literals)
                :id_(id),
                 graph_pattern_(graph_pattern),
                 new_vertexes_id_{new_vertex_id},
                 rhs_literals_(rhs_literals),
                 lhs_literals_(lhs_literals){
    return;
  }

  // construct method with multiply new vertexes
  ExpandTreeNode(int id,
                 const GraphPatternType& graph_pattern,
                 const std::vector<PatternVertexIDType>& new_vertexes_id,
                 const std::vector<DataGraphLiteralType>& rhs_literals,
                 const std::vector<DataGraphLiteralType>& lhs_literals)
                :id_(id),
                 graph_pattern_(graph_pattern),
                 new_vertexes_id_{new_vertexes_id},
                 rhs_literals_(rhs_literals),
                 lhs_literals_(lhs_literals) {
    return;
  }

  ExpandTreeNode(const ExpandTreeNode&) = default;
  ExpandTreeNode(ExpandTreeNode&&) = default;

  ExpandTreeNode& operator=(const ExpandTreeNode&) = default;	
  ExpandTreeNode& operator=(ExpandTreeNode&&) = default;	

  inline const int id() const{
    return this->id_;
  }

  inline GraphPatternType& pattern() {
    return this->graph_pattern_;
  }

  inline const GraphPatternType& pattern() const{
    return this->graph_pattern_;
  }
  
  inline const GraphPatternType& const_pattern() const{
    return this->graph_pattern_;
  }
  
  inline auto CountNewVertexes() const{
    return this->new_vertexes_id_.size();
  }

  inline const std::vector<PatternVertexIDType>& new_vertexes_id() const {
    return this->new_vertexes_id_;
  } 

  inline std::vector<PatternVertexIDType>& new_vertexes_id() {
    return this->new_vertexes_id_;
  } 

  inline typename std::vector<PatternVertexIDType>
                     ::const_iterator NewVertexesIdCBegin() const{
    assert(this->CountNewVertexes() == 0
        || this->CountNewVertexes() == 1
        || this->CountNewVertexes() == 2);
    return this->new_vertexes_id_.cbegin();
  }

  inline typename std::vector<PatternVertexIDType>
                     ::const_iterator NewVertexesIdCEnd() const{	
    assert(this->CountNewVertexes() == 0
        || this->CountNewVertexes() == 1
        || this->CountNewVertexes() == 2 // naive gcr implement
        );
    return this->new_vertexes_id_.cend();
  }

  inline typename std::vector<DataGraphLiteralType>
                     ::const_iterator RhsLiteralsBegin() const{
    return this->rhs_literals_.begin();
  }
  inline typename std::vector<DataGraphLiteralType>
                     ::const_iterator RhsLiteralsEnd() const{
    return this->rhs_literals_.end();
  }

  inline typename std::vector<DataGraphLiteralType>
                     ::iterator RhsLiteralsBegin(){
    return this->rhs_literals_.begin();
  }
  inline typename std::vector<DataGraphLiteralType>
                     ::iterator RhsLiteralsEnd(){
    return this->rhs_literals_.end();
  }

  inline typename std::vector<DataGraphLiteralType>
                     ::const_iterator LhsLiteralsBegin() const{
    return this->lhs_literals_.begin();
  }
  inline typename std::vector<DataGraphLiteralType>
                     ::const_iterator LhsLiteralsEnd() const{
    return this->lhs_literals_.end();
  }

  inline typename std::vector<DataGraphLiteralType>
                     ::iterator LhsLiteralsBegin(){
    return this->lhs_literals_.begin();
  }
  inline typename std::vector<DataGraphLiteralType>
                     ::iterator LhsLiteralsEnd(){
    return this->lhs_literals_.end();
  }

  inline const std::vector<DataGraphLiteralType>& const_rhs_literals() const{
    return this->rhs_literals_;
  }

  inline const std::vector<DataGraphLiteralType>& rhs_literals() const{
    return this->rhs_literals_;
  }

  inline std::vector<DataGraphLiteralType>& rhs_literals() {
    return this->rhs_literals_;
  }

  inline const std::vector<DataGraphLiteralType>& const_lhs_literals() const{
    return this->lhs_literals_;
  }

  inline const std::vector<DataGraphLiteralType>& lhs_literals() const{
    return this->lhs_literals_;
  }

  inline std::vector<DataGraphLiteralType>& lhs_literals() {
    return this->lhs_literals_;
  }

  inline typename std::vector<DataGraphLiteralType>::iterator 
     EraseRhsLiteral(typename std::vector<DataGraphLiteralType>::iterator literal_iterator){
    return this->rhs_literals_.erase(literal_iterator);
  }

  template<typename... Args>
  inline DataGraphLiteralType& AddRhsLiteral(const Args&... args){
    return this->rhs_literals_.emplace_back(args...);
  }

  inline void ClearRhsLiteral() {
    return this->rhs_literals_.clear();
  }

  template<typename... Args>
  inline DataGraphLiteralType& AddLhsLiteral(const Args&... args){
    return this->lhs_literals_.emplace_back(args...);
  }

 private:
  // the unique id in the same layer
  int id_;

  // the graph pattern
  GraphPatternType graph_pattern_;

  // the new added vertexes in this pattern
  // under current expanding strategy, it 
  // should only contain at most one element
  std::vector<PatternVertexIDType> new_vertexes_id_;

  // all possibel rhs_literals for this graph pattern
  // inherit from the parent node in the expanding tree
  //
  // according to the anti-monotonic if one literal cannot
  // satisfy the support bound before the pattern is expanded,
  // it cannot satisfy the support bound after it has been
  // expandation, therefore, in each expandation step,
  // it only needs to incrementally consider the new
  // rhs_literals introduced by the new vertexes
  std::vector<DataGraphLiteralType> rhs_literals_;

  // all lhs_literals for this graph pattern
  // inherit from the parent node in the expanding tree
  std::vector<DataGraphLiteralType> lhs_literals_;
};

// since the tree is expanded level-wise, there is no need to 
// preserve the entire tree, only need to preserve each layer
template<typename _GraphPatternType,
         typename     DataGraphType>
class ExpandTreeLevel {
 public:
  using GraphPatternType = _GraphPatternType;

  using DataGraphLiteralType = gar::LiteralInfo<_GraphPatternType,
                                                    DataGraphType>;

 private:
  using ExpandTreeNodeType 
      = ExpandTreeNode<GraphPatternType, 
                          DataGraphType>;

  using ExpandTreeLevelType
      = std::vector<ExpandTreeNodeType>;

 public:
  ExpandTreeLevel():level_(){
    return;
  }

  ExpandTreeLevel(const ExpandTreeLevel&) = default;
  ExpandTreeLevel(ExpandTreeLevel&&) = default;

  ExpandTreeLevel& operator=(const ExpandTreeLevel&) = default;	
  ExpandTreeLevel& operator=(ExpandTreeLevel&&) = default;	

  inline void Reserve(typename ExpandTreeLevelType::size_type size){
    this->level_.reserve(size);
    return;
  }

  // random access, for the open mp
  // it should be noticed that this idx might not align with the
  // id of the ExpandTreeNode 
  inline const ExpandTreeNodeType& node(
      typename ExpandTreeLevelType::size_type idx) const{
    assert(this->level_.size() > idx);
    return this->level_[idx];
  }
  inline ExpandTreeNodeType& node(
      typename ExpandTreeLevelType::size_type idx){
    assert(this->level_.size() > idx);
    return this->level_[idx];
  }

  inline const ExpandTreeNodeType& back() const{
    return this->level_.back();
  }
  
  inline ExpandTreeNodeType& back(){
    return this->level_.back();
  }

  inline void pop_back() {
    this->level_.pop_back();
    return;
  }

  inline typename ExpandTreeLevelType::size_type size() const {
    return this->level_.size();
  }

  template<typename... Args>
  inline ExpandTreeNodeType& AddExpandTreeNode(const Args&... args){
    return this->level_.emplace_back(args...);
  }

  inline void MoveToEnd(ExpandTreeLevel& expand_level){
    this->level_.insert(this->level_.end(),
         std::make_move_iterator(expand_level.level_.begin()),
         std::make_move_iterator(expand_level.level_.end()));
    return;
  }

  inline typename ExpandTreeLevelType::const_iterator begin() const{
    return this->level_.begin();
  }
  inline typename ExpandTreeLevelType::const_iterator end() const{
    return this->level_.end();
  }

  inline typename ExpandTreeLevelType::iterator begin(){
    return this->level_.begin();
  }
  inline typename ExpandTreeLevelType::iterator end(){
    return this->level_.end();
  }

  inline void clear() noexcept {
    this->level_.clear();
    return;
  }

  inline bool empty() const {
    return this->level_.empty();
  }

  inline void Swap(ExpandTreeLevel& expand_level){
    this->level_.swap(expand_level.level_);
    return;
  }

  // set the id of nodes same as the index of it
  inline void ReorderNodeId(){
    for (int i = 0; i < this->level_.size(); i++){
      this->level_[i].set_id(i);
    }
    return;
  }

 private:
  ExpandTreeLevelType level_;
};

}  // namespace _gar_discover

}  // namespace grape

#endif // _EXPAND_TREE_H