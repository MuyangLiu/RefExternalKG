#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_RESTRICTION_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_RESTRICTION_H_

#include "util/log.h"

namespace grape {

namespace _gar_discover {

template<typename GraphPatternType,
         typename    DataGraphType>
class Restriction {
 private:
  using LiteralStandAloneInfoType
  =gar::LiteralStandAloneInfo<GraphPatternType,
                                 DataGraphType>;

  using DataGraphVertexIDType 
      = typename GUNDAM::VertexID<DataGraphType>::type;

  using DataGraphVertexHandleType 
      = typename GUNDAM::VertexHandle<DataGraphType>::type;

 public:
  Restriction() : /* #################### *
                   * ##  restrictions  ## *
                   * #################### */
                  pattern_is_link_(false),
                  pattern_is_star_(false),
                  pattern_is_tree_(false),
                  edge_literal_only_between_2_hop_connected_vertexes_(false),
              variable_literal_only_between_2_hop_connected_vertexes_(false),
              variable_literal_only_between_connected_vertexes_(false), 
                       literals_connected_(false),
                        discretized_graph_(false),
                       central_to_central_(false),
                  /* ################## *
                   * ##  rule types  ## *
                   * ################## */
                  horn_rule_(false),
                        gar_(false),
                        gcr_(false),
                     pattern_vertex_limit_(0),
                           diameter_limit_(0),
                           path_num_limit_(0), // for gcr
                        path_length_limit_(0), // for gcr
                       bidirectional_path_(true), // for gcr
                consider_constant_literal_(true),
                consider_variable_literal_(true),
                    consider_edge_literal_(true),
                      consider_ml_literal_(true),
         consider_constant_literal_in_lhs_(true),
         consider_variable_literal_in_lhs_(true),
               consider_ml_literal_in_lhs_(true),
                      constant_freq_bound_(0.0),
               specified_confidence_bound_(false),
                         confidence_bound_(0.0) {
    return;
  }

  inline bool AddRestriction(const std::string& restriction) {
    if (restriction == "horn_rule") {
      util::Error("restriction: horn_rule have been moved to SpecifyRuleType");
      return false;
      // if (this->pattern_is_link_) {
      //   util::Error("conflict with restriction: pattern_link_only");
      //   return false;
      // }
      // if (this->pattern_is_link_) {
      //   util::Error("conflict with restriction: pattern_without_loop");
      //   return false;
      // }
      // if (this->edge_literal_only_between_2_hop_connected_vertexes_) {
      //   util::Error("conflict with restriction: edge_literal_only_between_2_hop_connected_vertexes");
      //   return false;
      // }
      // if (this->variable_literal_only_between_connected_vertexes_) {
      //   util::Error("conflict with restriction: variable_literal_only_between_connected_vertexes");
      //   return false;
      // }
      // if (this->literals_connected_) {
      //   util::Error("conflict with restriction: literals_connected");
      //   return false;
      // }
      // this->horn_rule_ = true;
      // return true;
    }
    if (restriction == "pattern_is_star") {
      util::Info("with restriction: pattern_is_star");
      this->pattern_is_star_ = true;
      return true;
    }
    if (restriction == "pattern_without_loop"
     || restriction == "pattern_is_tree") {
      util::Info("with restriction: pattern_without_loop ");
      // if (this->horn_rule_) {
      //   util::Error("conflict with restriction: horn_rule");
      //   return false;
      // }
      if (this->pattern_is_tree_) {
        util::Error("conflict with restriction: pattern_is_link");
        return false;
      }
      this->pattern_is_tree_ = true;
      return true;
    }
    if (restriction == "edge_literal_only_between_2_hop_connected_vertexes") {
      util::Info("with restriction: edge_literal_only_between_2_hop_connected_vertexes");
      // if (this->horn_rule_) {
      //   util::Error("conflict with restriction: horn_rule");
      //   return false;
      // }
      this->edge_literal_only_between_2_hop_connected_vertexes_ = true;
      return true;
    }
    if (restriction == "variable_literal_only_between_2_hop_connected_vertexes") {
      util::Info("with restriction: variable_literal_only_between_2_hop_connected_vertexes");
      // if (this->horn_rule_) {
      //   util::Error("conflict with restriction: horn_rule");
      //   return false;
      // }
      if (this->variable_literal_only_between_connected_vertexes_) {
        util::Error("conflict with restriction: variable_literal_only_between_connected_vertexes");
        return false;
      }
      this->variable_literal_only_between_2_hop_connected_vertexes_ = true;
      return true;
    }
    if (restriction == "variable_literal_only_between_connected_vertexes") {
      util::Info("with restriction: variable_literal_only_between_connected_vertexes");
      // if (this->horn_rule_) {
      //   util::Error("conflict with restriction: horn_rule");
      //   return false;
      // }
      if (this->variable_literal_only_between_2_hop_connected_vertexes_) {
        util::Error("conflict with restriction: variable_literal_only_between_2_hop_connected_vertexes");
        return false;
      }
      this->variable_literal_only_between_connected_vertexes_ = true;
      return true;
    }
    if (restriction == "literals_connected") {
      util::Info("with restriction: literals_connected");
      // if (this->horn_rule_) {
      //   util::Error("conflict with restriction: horn_rule");
      //   return false;
      // }
      this->literals_connected_ = true;
      return true;
    }
    if (restriction == "discretized_graph") {
      if (this->discretized_graph_) {
        util::Error("duplciated restriction: discretized_graph ");
        return false;
      }
      util::Info("with restriction: discretized_graph ");
      this->discretized_graph_ = true;
      return true;
    }
    if (restriction == "central_to_central") {
      if (this->central_to_central_) {
        util::Error("duplciated restriction: central_to_central ");
        return false;
      }
      util::Info("with restriction: central_to_central ");
      this->central_to_central_ = true;
      return true;
    }
    util::Error("unknown restriction: " + restriction);
    return false;
  }

  inline const bool pattern_is_tree() const {
    return this->pattern_is_tree_;
  }

  inline const bool pattern_is_link() const {
    return this->pattern_is_link_;
  }

  inline const bool pattern_is_star() const {
    return this->pattern_is_star_;
  }
  
  inline const bool edge_literal_only_between_2_hop_connected_vertexes() const {
    return this->edge_literal_only_between_2_hop_connected_vertexes_;
  }

  inline const bool variable_literal_only_between_2_hop_connected_vertexes() const {
    return this->variable_literal_only_between_2_hop_connected_vertexes_;
  }

  inline const bool variable_literal_only_between_connected_vertexes() const {
    return this->variable_literal_only_between_connected_vertexes_;
  }

  inline const bool literals_connected() const {
    return this->literals_connected_;
  }

  inline bool set_pattern_vertex_limit(size_t pattern_vertex_limit) {
    if (pattern_vertex_limit <= 0) {
      util::Error("illegal pattern_vertex_limit: "
          + std::to_string(pattern_vertex_limit));
      return false;
    }
    if (this->pattern_vertex_limit_ != 0) {
      util::Error("pattern_vertex_limit has been set");
      return false;
    }
    this->pattern_vertex_limit_ = pattern_vertex_limit;
    return true;
  }

  inline size_t pattern_vertex_limit() const {
    assert(this->pattern_vertex_limit_ > 0);
    return this->pattern_vertex_limit_;
  }

  inline bool has_pattern_vertex_limit() const {
    return this->pattern_vertex_limit_ > 0;
  }

  inline bool set_diameter_limit(size_t diameter_limit) {
    if (diameter_limit <= 0) {
      util::Error("illegal diameter_limit: "
          + std::to_string(diameter_limit));
      return false;
    }
    if (this->diameter_limit_ != 0) {
      util::Error("diameter_limit has been set");
      return false;
    }
    this->diameter_limit_ = diameter_limit;
    return true;
  }

  inline size_t diameter_limit() const {
    assert(this->diameter_limit_ > 0);
    return this->diameter_limit_;
  }

  inline bool has_diameter_limit() const {
    return this->diameter_limit_ > 0;
  }

  inline void set_consider_no_literal() {
    this->consider_constant_literal_ = false;
    this->consider_variable_literal_ = false;
    this->    consider_edge_literal_ = false;
    this->      consider_ml_literal_ = false;
    return;
  }

  inline void set_consider_no_literal_in_lhs() {
    this->consider_constant_literal_in_lhs_ = false;
    this->consider_variable_literal_in_lhs_ = false;
    this->      consider_ml_literal_in_lhs_ = false;
    return;
  }

  inline bool discretized_graph() const {
    return this->discretized_graph_;
  }

  inline bool central_to_central() const {
    return this->central_to_central_;
  }

  inline bool ConsiderLiteral(const std::string& literal_type) {
    if (literal_type == "constant_literal") {
      if (this->consider_constant_literal_) {
        util::Error("duplciated literal type: constant_literal ");
        return false;
      }
      util::Info("with literal type: constant_literal ");
      this->consider_constant_literal_ = true;
      return true;
    }
    if (literal_type == "variable_literal") {
      if (this->consider_variable_literal_) {
        util::Error("duplciated literal type: variable_literal ");
        return false;
      }
      util::Info("with literal type: variable_literal");
      this->consider_variable_literal_ = true;
      return true;
    }
    if (literal_type == "edge_literal") {
      if (this->consider_edge_literal_) {
        util::Error("duplciated literal type: edge_literal ");
        return false;
      }
      util::Info("with literal type: edge_literal");
      this->consider_edge_literal_ = true;
      return true;
    }
    if (literal_type == "ml_literal") {
      if (this->consider_ml_literal_) {
        util::Error("duplciated literal type: ml_literal ");
        return false;
      }
      util::Info("with literal type: ml_literal");
      this->consider_ml_literal_ = true;
      return true;
    }
    util::Error("unknown literal type: " + literal_type);
    return false;
  }

  inline bool SetPatternType(const std::string& pattern_type) {
    if (pattern_type == "star") {
      this->pattern_is_star_ = true;
      return true;
    }
    if (pattern_type == "tree") {
      this->pattern_is_tree_ = true;
      return true;
    }
    if (pattern_type == "link") {
      this->pattern_is_link_ = true;
      return true;
    }
    util::Error("unknown pattern type: " 
                       + pattern_type);
    return false;
  }

  inline bool ConsiderLiteralInLhs(const std::string& literal_type) {
    if (literal_type == "constant_literal") {
      if (this->consider_constant_literal_in_lhs_) {
        util::Error("duplciated literal type: constant_literal ");
        return false;
      }
      if (!this->consider_constant_literal_) {
        util::Error("does not consider constant_literal at all, but add it to Lhs");
        return false;
      }
      util::Info("with literal type in lhs: constant_literal ");
      this->consider_constant_literal_in_lhs_ = true;
      return true;
    }
    if (literal_type == "variable_literal") {
      if (this->consider_variable_literal_in_lhs_) {
        util::Error("duplciated literal type: variable_literal ");
        return false;
      }
      if (!this->consider_variable_literal_) {
        util::Error("does not consider variable_literal at all, but add it to Lhs");
        return false;
      }
      util::Info("with literal type in lhs: variable_literal");
      this->consider_variable_literal_in_lhs_ = true;
      return true;
    }
    if (literal_type == "edge_literal") {
      util::Error("does not support edge_literal in lhs!");
      return false;
    }
    if (literal_type == "ml_literal") {
      if (this->consider_ml_literal_in_lhs_) {
        util::Error("duplciated literal type: ml_literal ");
        return false;
      }
      if (!this->consider_ml_literal_) {
        util::Error("does not consider ml_literal at all, but add it to Lhs");
        return false;
      }
      util::Info("with literal type in lhs: ml_literal");
      this->consider_ml_literal_in_lhs_ = true;
      return true;
    }
    util::Error("unknown literal type in lhs: " + literal_type);
    return false;
  }

  inline bool LiteralTypeConsideredInLhs(const gar::LiteralType& literal_type) const {
    switch(literal_type) {
     case gar::LiteralType::kConstantLiteral:
      return this->consider_constant_literal_in_lhs_;
     case gar::LiteralType::kVariableLiteral:
      return this->consider_variable_literal_in_lhs_;
     case gar::LiteralType::kMlLiteral:
      return this->consider_ml_literal_in_lhs_;
     case gar::LiteralType::kEdgeLiteral:
      // does not consider edge literal in lhs
      return false;
     default:
      assert(false);
    }
    return false;
  }

  inline bool SpecifyRuleType(const std::string& rule_type) {
    if (rule_type == "horn_rule") {
      util::Info("rule type: horn_rule");
      if (this->horn_rule_
       || this->gar_
       || this->gcr_) {
        util::Error("duplicated rule type specification");
        return false;
      }
      this->horn_rule_ = true;
      this->pattern_is_link_ = true;
      return true;
    }
    if (rule_type == "gar") {
      util::Info("rule type: gar");
      if (this->horn_rule_
       || this->gar_
       || this->gcr_) {
        util::Error("duplicated rule type specification");
        return false;
      }
      this->gar_ = true;
      return true;
    }
    if (rule_type == "gcr") {
      util::Info("rule type: gcr");
      if (this->horn_rule_
       || this->gar_
       || this->gcr_) {
        util::Error("duplicated rule type specification");
        return false;
      }
      this->gcr_ = true;
      if (!this->pattern_is_link() 
       && !this->pattern_is_tree() 
       && !this->pattern_is_star() ) {
        this->pattern_is_star_ = true;
      }
      return true;
    }
    util::Error("unknown rule type: " + rule_type);
    return false;
  }

  /* ################## *
   * ##  rule types  ## *
   * ################## */
  inline const bool horn_rule() const {
    return this->horn_rule_;
  }
  inline const bool gcr() const {
    return this->gcr_;
  }
  inline const bool gar() const {
    return this->gar_;
  }

  inline bool consider_constant_literal() const {
    return this->consider_constant_literal_;
  }

  inline bool consider_variable_literal() const {
    return this->consider_variable_literal_;
  }

  inline bool consider_edge_literal() const {
    return this->consider_edge_literal_;
  }

  inline bool consider_ml_literal() const {
    return this->consider_ml_literal_;
  }

  inline void set_constant_freq_bound(float constant_freq_bound) {
    this->constant_freq_bound_
        = constant_freq_bound;
    return;
  }

  inline float constant_freq_bound() const {
    return this->constant_freq_bound_;
  }

  inline void set_confidence_bound(float confidence_bound) {
    this->specified_confidence_bound_ = true;
    this->confidence_bound_
        = confidence_bound;
    return;
  }

  inline void set_path_num_limit(int path_num_limit) {
    assert(!this->specified_path_num_limit());
    assert(path_num_limit >= 1);
    this-> path_num_limit_ = path_num_limit;
    assert( this->specified_path_num_limit());
    return;
  }

  inline bool specified_path_num_limit() const {
    return this->path_num_limit_ >= 1;
  }

  inline size_t path_num_limit() const {
    assert(this->specified_path_num_limit());
    return this->path_num_limit_;
  }

  inline void set_path_length_limit(int path_length_limit) {
    assert(!this->specified_path_length_limit());
    assert(path_length_limit >= 1);
    this-> path_length_limit_ = path_length_limit;
    assert( this->specified_path_length_limit());
    return;
  }

  inline bool specified_path_length_limit() const {
    return this->path_length_limit_ >= 1;
  }

  inline size_t path_length_limit() const {
    assert(this->specified_path_length_limit());
    return this->path_length_limit_;
  }

  inline bool specified_confidence_bound() const {
    return this->specified_confidence_bound_;
  }

  inline float confidence_bound() const {
    return this->confidence_bound_;
  }

  inline void set_bidirectional_path(bool kDirectedPath) {
    this->bidirectional_path_ = kDirectedPath;
    return;
  }

  inline bool bidirectional_path() const {
    return this->bidirectional_path_;
  }

  inline bool AddSpecifiedRhsLiteral(
      const LiteralStandAloneInfoType& literal_stand_alone_info,
      const std::vector<
            std::vector<DataGraphVertexHandleType>>& literal_supp_set = {}) {
    std::vector<
    std::vector<DataGraphVertexHandleType>> candidate_set_for_vertex_in_rhs_literal;
    if (!literal_supp_set.empty()) {
      candidate_set_for_vertex_in_rhs_literal.resize(literal_supp_set.begin()->size());
      for (const auto& literal_supp : literal_supp_set) {
        for (size_t literal_supp_idx = 0;
                    literal_supp_idx < literal_supp.size();
                    literal_supp_idx++) {
          assert(literal_supp_idx < candidate_set_for_vertex_in_rhs_literal.size());
          candidate_set_for_vertex_in_rhs_literal[literal_supp_idx]
                   .emplace_back(literal_supp[literal_supp_idx]);
        }
      }
      for (auto& candidate_set 
               : candidate_set_for_vertex_in_rhs_literal) {
        std::sort(candidate_set.begin(),
                  candidate_set.end());
        candidate_set.erase(std::unique(candidate_set.begin(), 
                                        candidate_set.end()),
                                        candidate_set.end());
      }
    }
    return this->specified_rhs_literal_set_
                .emplace(literal_stand_alone_info,
               std::pair(literal_supp_set,
                         candidate_set_for_vertex_in_rhs_literal)).second;
  }

  inline void AddSpecifiedPath(const GraphPatternType& path) {
    this->specified_path_set_
         .emplace_back(path);
    return;
  }

  inline const auto& specified_path_set() const {
    return this->specified_path_set_;
  }

  inline const auto& specified_rhs_literal_set() const {
    return this->specified_rhs_literal_set_;
  }

  // specified this rhs literal and its candidate set is not empty
  inline bool specified_support_set_for_rhs_literal(
           const LiteralStandAloneInfoType& rhs_literal_stand_alone_info) const {
       auto specified_rhs_literal_set_it
    = this->specified_rhs_literal_set_.find(rhs_literal_stand_alone_info);
    if (specified_rhs_literal_set_it
      == this->specified_rhs_literal_set_.end())  {
      // does not specified this rhs literal
      return false;
    }
    return !specified_rhs_literal_set_it->second.first.empty();
  }

  inline const auto& support_set_for_rhs_literal(
                      const LiteralStandAloneInfoType& rhs_literal_stand_alone_info) const {
    assert(this->specified_support_set_for_rhs_literal(rhs_literal_stand_alone_info));
    return this->specified_rhs_literal_set_
                .find(rhs_literal_stand_alone_info)->second.first;
  }

  inline const auto& candidate_set_for_vertex_in_rhs_literal(
                      const LiteralStandAloneInfoType& rhs_literal_stand_alone_info) const {
    assert(this->specified_support_set_for_rhs_literal(rhs_literal_stand_alone_info));
    return this->specified_rhs_literal_set_
                .find(rhs_literal_stand_alone_info)->second.second;
  }

 private:
  // for restirctions
  bool pattern_is_link_,
       pattern_is_star_,
       pattern_is_tree_,
       edge_literal_only_between_2_hop_connected_vertexes_,
   variable_literal_only_between_2_hop_connected_vertexes_,
         variable_literal_only_between_connected_vertexes_,
         literals_connected_,
         // for gcr
          discretized_graph_,
         central_to_central_;

  // for rule_type
  bool horn_rule_,
             gar_,
             gcr_;

  size_t pattern_vertex_limit_,
               diameter_limit_;

  // to mark the type of literals to be consider
  bool consider_constant_literal_,
       consider_variable_literal_,
           consider_edge_literal_,
             consider_ml_literal_;

  // to mark the type of literals to be consider in lhs
  bool consider_constant_literal_in_lhs_,
       consider_variable_literal_in_lhs_,
             consider_ml_literal_in_lhs_;

  // consider only the frequent value
  float constant_freq_bound_;

  // for confidence bound
  bool specified_confidence_bound_;
  float confidence_bound_;

  std::map<LiteralStandAloneInfoType,
  std::pair<
  std::vector<
  std::vector<DataGraphVertexHandleType>>, // support_set_for_rhs_literal
  std::vector<
  std::vector<DataGraphVertexHandleType>>>> specified_rhs_literal_set_;

  // for gcr
  std::vector<GraphPatternType> specified_path_set_;

  int path_num_limit_,
   path_length_limit_;

  bool bidirectional_path_;
};

}  // namespace _gar_discover

}  // namespace grape

#endif // EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_RESTRICTION_H_