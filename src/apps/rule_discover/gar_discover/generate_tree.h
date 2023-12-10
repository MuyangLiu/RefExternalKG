#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GENERATE_TREE_H
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GENERATE_TREE_H

#include "gundam/graph_type/large_graph2.h"
#include "gar/literal.h"
#include "gar/gar.h"

#include "rule_discover/gar_discover/expand_tree.h"

namespace grape {

namespace _gar_discover {

using SubstructureLevelType = int;

//using LiteralAttributeKeyType = bool;

using LiteralAttributeKeyType = int;

using LiteralTreeType
    = GUNDAM::LargeGraph2<int, SubstructureLevelType, LiteralAttributeKeyType,
                          int, int, int>;

using LiteralTreePtr = typename GUNDAM::VertexHandle<LiteralTreeType>::type;

using LiteralTreeDefaultEdgeLabelType
  = typename _gar_discover::LiteralTreeType 
                                 ::EdgeType 
                                ::LabelType;

LiteralTreeDefaultEdgeLabelType kLiteralTreeDefaultEdgeLabel = 0;

//static constexpr LiteralAttributeKeyType kLiteralKey = true;
//static constexpr LiteralAttributeKeyType    kConfKey = false;

static constexpr LiteralAttributeKeyType kLiteralKey = 1;
static constexpr LiteralAttributeKeyType    kConfKey = 0;
static constexpr LiteralAttributeKeyType    kSuppKey = 2;
static constexpr LiteralAttributeKeyType kNewSuppKey = 3;
static constexpr LiteralAttributeKeyType kMatchFileName = 4;
static constexpr LiteralAttributeKeyType kOriginalGarID = 5;
static constexpr LiteralAttributeKeyType kProbability = 6;

template<typename GraphPatternType,
         typename  LiteralTreeType>
class GenerateTreeNode;

template<typename GraphPatternType,
         typename    DataGraphType>
class GenerateTreeLevel;

template<typename GraphPatternType,
         typename    DataGraphType>
class GenerateTreeNode{
 private:
  using LiteralInfoType
          = gar::LiteralInfo<GraphPatternType,
                               DataGraphType>;
                               
  using LiteralStandAloneInfoType
          = gar::LiteralStandAloneInfo<GraphPatternType,
                                          DataGraphType>;
 public:
  GenerateTreeNode(int id, 
                   const GraphPatternType& graph_pattern)
                  :id_(id),
                   graph_pattern_(graph_pattern){
    return;
  }
  
  // forbidden copy constructure method
  GenerateTreeNode(const GenerateTreeNode&) = delete;
  GenerateTreeNode& operator=(const GenerateTreeNode&) = delete;	
  
  // only allow the move constructure
  GenerateTreeNode(GenerateTreeNode&&) = default;
  GenerateTreeNode& operator=(GenerateTreeNode&&) = default;	

  inline int id() const{
    return this->id_;
  }

  inline GraphPatternType& pattern(){
    return this->graph_pattern_;
  }

  inline const GraphPatternType& const_pattern() const{
    return this->graph_pattern_;
  }

  inline void AddLiteralTree(const LiteralInfoType& literal,
                             const float            literal_conf) {
    auto& literal_tree_ref = this->literal_trees_.emplace_back();
    // <vertex_ptr, bool>
    auto add_vertex_ret = literal_tree_ref.AddVertex(0, 1);
    // should added successfully
    assert(add_vertex_ret.first);
    assert(add_vertex_ret.second);
    // kLiteralKey
    auto add_literal_attr_ret
       = add_vertex_ret.first->AddAttribute(kLiteralKey, literal);
    assert(!add_literal_attr_ret.first.IsNull());
    assert( add_literal_attr_ret.second);
    // kConfKey
    auto add_conf_attr_ret
       = add_vertex_ret.first->AddAttribute(kConfKey, literal_conf);
    assert(!add_conf_attr_ret.first.IsNull());
    assert(add_conf_attr_ret.second);
    return;
  }

  inline void AddLiteralTreeWithSupp(const LiteralInfoType& literal,
                             const float            literal_conf,
                             const int              literal_supp,
                             const int              literal_gar_id,
                             const std::string&     literal_pivot_file_name,
                             const double           probability) {
    auto& literal_tree_ref = this->literal_trees_.emplace_back();
    // <vertex_ptr, bool>
    auto add_vertex_ret = literal_tree_ref.AddVertex(0, 1);
    // should added successfully
    assert(add_vertex_ret.first);
    assert(add_vertex_ret.second);
    // kLiteralKey
    auto add_literal_attr_ret
       = add_vertex_ret.first->AddAttribute(kLiteralKey, literal);
    assert(!add_literal_attr_ret.first.IsNull());
    assert( add_literal_attr_ret.second);
    // kConfKey
    auto add_conf_attr_ret
       = add_vertex_ret.first->AddAttribute(kConfKey, literal_conf);
    assert(!add_conf_attr_ret.first.IsNull());
    assert(add_conf_attr_ret.second);

    auto add_supp_attr_ret
       = add_vertex_ret.first->AddAttribute(kSuppKey, literal_supp);
    assert(!add_supp_attr_ret.first.IsNull());
    assert(add_supp_attr_ret.second);
    auto add_id_attr_ret
       = add_vertex_ret.first->AddAttribute(kOriginalGarID, literal_gar_id);
    assert(!add_id_attr_ret.first.IsNull());
    assert(add_id_attr_ret.second);
    auto pivot_file_ret
       = add_vertex_ret.first->AddAttribute(kMatchFileName, literal_pivot_file_name);
    assert(!pivot_file_ret.first.IsNull());
    assert(pivot_file_ret.second);
    auto probability_ret
       = add_vertex_ret.first->AddAttribute(kProbability, probability);
    assert(!probability_ret.first.IsNull());
    assert(probability_ret.second);
    return;
  }


  inline void InitLiteralTrees(const std::vector<LiteralInfoType>& literals_set,
                               const std::vector<float>&           literals_conf_set) {
    assert(literals_set.size() == literals_conf_set.size());
    this->literal_trees_.reserve( literals_set.size() );
    this->literals_set_num_ = literals_set.size();
    // for all possible literal, build a tree with the root of it
    for (size_t literal_idx = 0;
                literal_idx < literals_set.size();
                literal_idx++) {
      const auto& literal      = literals_set     [literal_idx];
      const auto& literal_conf = literals_conf_set[literal_idx];
      // construct a new literal tree
      auto& literal_tree_ref = this->literal_trees_.emplace_back();
      // initialize the literal trees
      // <vertex_ptr, bool>
      auto add_vertex_ret = literal_tree_ref.AddVertex(0, 1);
      // should added successfully
      assert(add_vertex_ret.first);
      assert(add_vertex_ret.second);
      // kLiteralKey
      auto add_literal_attr_ret
         = add_vertex_ret.first->AddAttribute(kLiteralKey, literal);
      assert(!add_literal_attr_ret.first.IsNull());
      assert(add_literal_attr_ret.second);
      // kConfKey
      auto add_conf_attr_ret
         = add_vertex_ret.first->AddAttribute(kConfKey, literal_conf);
      assert(!add_conf_attr_ret.first.IsNull());
      assert(add_conf_attr_ret.second);
    }
    return;
  }


  inline void InitLiteralTreesWithSupp(const std::vector<LiteralInfoType>& literals_set,
                               const std::vector<float>&           literals_conf_set,
                               const std::vector<int>&             literals_supp_set) {
    assert(literals_set.size() == literals_conf_set.size());
    this->literal_trees_.reserve( literals_set.size() );
    this->literals_set_num_ = literals_set.size();
    // for all possible literal, build a tree with the root of it
    for (size_t literal_idx = 0;
                literal_idx < literals_set.size();
                literal_idx++) {
      const auto& literal      = literals_set     [literal_idx];
      const auto& literal_conf = literals_conf_set[literal_idx];
      const auto& literal_supp = literals_supp_set[literal_idx];

      // construct a new literal tree
      auto& literal_tree_ref = this->literal_trees_.emplace_back();
      // initialize the literal trees
      // <vertex_ptr, bool>
      auto add_vertex_ret = literal_tree_ref.AddVertex(0, 1);
      // should added successfully
      assert(add_vertex_ret.first);
      assert(add_vertex_ret.second);
      // kLiteralKey
      auto add_literal_attr_ret
         = add_vertex_ret.first->AddAttribute(kLiteralKey, literal);
      assert(!add_literal_attr_ret.first.IsNull());
      assert(add_literal_attr_ret.second);
      // kConfKey
      auto add_conf_attr_ret
         = add_vertex_ret.first->AddAttribute(kConfKey, literal_conf);
      assert(!add_conf_attr_ret.first.IsNull());
      assert(add_conf_attr_ret.second);

      // kSuppKey
      auto add_supp_attr_ret
         = add_vertex_ret.first->AddAttribute(kSuppKey, literal_supp);
      assert(!add_supp_attr_ret.first.IsNull());
      assert(add_supp_attr_ret.second);
    }
    return;
  }

  inline void InitLiteralTreesWithSuppAndName(const std::vector<LiteralInfoType>& literals_set,
                               const std::vector<float>&           literals_conf_set,
                               const std::vector<int>&             literals_supp_set,
                               const std::vector<std::string>&     literals_name_set,
                               const std::vector<double>& literals_probability_set) {
    assert(literals_set.size() == literals_conf_set.size());
    this->literal_trees_.reserve( literals_set.size() );
    this->literals_set_num_ = literals_set.size();
    // for all possible literal, build a tree with the root of it
    for (size_t literal_idx = 0;
                literal_idx < literals_set.size();
                literal_idx++) {
      const auto& literal      = literals_set     [literal_idx];
      const auto& literal_conf = literals_conf_set[literal_idx];
      const auto& literal_supp = literals_supp_set[literal_idx];
      const auto& literal_name = literals_name_set[literal_idx];
      const auto& probability  = literals_probability_set[literal_idx];

      // construct a new literal tree
      auto& literal_tree_ref = this->literal_trees_.emplace_back();
      // initialize the literal trees
      // <vertex_ptr, bool>
      auto add_vertex_ret = literal_tree_ref.AddVertex(0, 1);
      // should added successfully
      assert(add_vertex_ret.first);
      assert(add_vertex_ret.second);
      // kLiteralKey
      auto add_literal_attr_ret
         = add_vertex_ret.first->AddAttribute(kLiteralKey, literal);
      assert(!add_literal_attr_ret.first.IsNull());
      assert(add_literal_attr_ret.second);
      // kConfKey
      auto add_conf_attr_ret
         = add_vertex_ret.first->AddAttribute(kConfKey, literal_conf);
      assert(!add_conf_attr_ret.first.IsNull());
      assert(add_conf_attr_ret.second);

      // kSuppKey
      auto add_supp_attr_ret
         = add_vertex_ret.first->AddAttribute(kSuppKey, literal_supp);
      assert(!add_supp_attr_ret.first.IsNull());
      assert(add_supp_attr_ret.second);

      auto add_name_attr_ret
         = add_vertex_ret.first->AddAttribute(kMatchFileName, literal_name);
      assert(!add_name_attr_ret.first.IsNull());
      assert(add_name_attr_ret.second);

      auto add_probability_ret
         = add_vertex_ret.first->AddAttribute(kProbability, probability);
      assert(!add_probability_ret.first.IsNull());
      assert(add_probability_ret.second);
    }
    return;
  }


  inline typename std::vector<LiteralTreeType>::iterator ExpandableLiteralTreesBegin() {
    return this->literal_trees_.begin();
  }

  inline typename std::vector<LiteralTreeType>::iterator ExpandableLiteralTreesEnd() {
    return this->literal_trees_.begin()
         + this->literals_set_num_;
  }

  inline typename std::vector<LiteralTreeType>::const_iterator ExpandableLiteralTreesBegin() const {
    return this->literal_trees_.begin();
  }

  inline typename std::vector<LiteralTreeType>::const_iterator ExpandableLiteralTreesEnd() const {
    return this->literal_trees_.begin()
         + this->literals_set_num_;
  }

  inline typename std::vector<LiteralTreeType>::iterator LiteralTreesBegin() {
    return this->literal_trees_.begin();
  }

  inline typename std::vector<LiteralTreeType>::iterator LiteralTreesEnd() {
    return this->literal_trees_.end();
  }

  inline typename std::vector<LiteralTreeType>::const_iterator LiteralTreesBegin() const {
    return this->literal_trees_.begin();
  }

  inline typename std::vector<LiteralTreeType>::const_iterator LiteralTreesEnd() const {
    return this->literal_trees_.end();
  }

  inline typename std::vector<LiteralTreeType>::iterator LiteralTreesBackIter() {
    if (this->literal_trees_.size() == 0) {
      return this->literal_trees_.begin();
    }
    return this->literal_trees_.end() - 1; 
  }

  inline typename std::vector<LiteralTreeType>::size_type LiteralTreesSize() const {
    return this->literal_trees_.size();
  }

  inline LiteralTreeType LiteralTreesBack() {
    return this->literal_trees_.back();
  }

 private:
  int id_;

  GraphPatternType graph_pattern_;

  std::vector<LiteralTreeType> literal_trees_;

  size_t literals_set_num_;
};

template<typename GraphPatternType,
         typename    DataGraphType>
class GenerateTreeLevel{
 private:
  using GenerateTreeNodeType
      = GenerateTreeNode<GraphPatternType,
                            DataGraphType>;

  using GenerateTreeLevelType
     = std::vector<GenerateTreeNodeType>;

 public:
  GenerateTreeLevel():level_(){
    return;
  }

  // forbidden copy constructure method  
  GenerateTreeLevel(const GenerateTreeLevel&) = delete;
  GenerateTreeLevel& operator=(const GenerateTreeLevel&) = delete;	

  // only allow the move constructure
  GenerateTreeLevel(GenerateTreeLevel&&) = default;
  GenerateTreeLevel& operator=(GenerateTreeLevel&&) = default;	
  
  // random access, for the open mp
  // it should be noticed that this idx might not align with the
  // id of the ExpandTreeNode 
  inline const GenerateTreeNodeType& node(
      typename GenerateTreeLevelType::size_type idx) const{
    assert(this->level_.size() > idx);
    return this->level_[idx];
  }
  inline GenerateTreeNodeType& node(
      typename GenerateTreeLevelType::size_type idx){
    assert(this->level_.size() > idx);
    return this->level_[idx];
  }

  template<typename... Args>
  inline GenerateTreeNodeType& AddGenerateTreeNode(Args&&... args){
    return this->level_.emplace_back(std::forward<Args...>(args...));
  }

  inline void MoveToEnd(GenerateTreeLevel& generate_level){
    this->level_.insert(this->level_.end(),
         std::make_move_iterator(generate_level.level_.begin()),
         std::make_move_iterator(generate_level.level_.end()));
    return;
  }

  inline typename GenerateTreeLevelType::const_iterator begin() const{
    return this->level_.begin();
  }
  inline typename GenerateTreeLevelType::const_iterator end() const{
    return this->level_.end();
  }

  inline typename GenerateTreeLevelType::iterator begin(){
    return this->level_.begin();
  }
  inline typename GenerateTreeLevelType::iterator end(){
    return this->level_.end();
  }

  inline void clear() noexcept {
    this->level_.clear();
    return;
  }

  inline typename GenerateTreeLevelType::size_type size() const {
    return this->level_.size();
  }

  inline void Reserve(typename GenerateTreeLevelType::size_type size){
    this->level_.reserve(size);
    return;
  }

 private:
  GenerateTreeLevelType level_;
};


}  // namespace _gar_discover

}  // namespace grape

#endif // EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GENERATE_TREE_H
