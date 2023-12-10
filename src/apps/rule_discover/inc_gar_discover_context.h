#ifndef EXAMPLES_ANALYTICAL_APPS_INC_GAR_DISCOVER_GAR_DISCOVER_CONTEXT_H_
#define EXAMPLES_ANALYTICAL_APPS_INC_GAR_DISCOVER_GAR_DISCOVER_CONTEXT_H_

#include <grape/app/context_base.h>
#include <grape/grape.h>

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "gar/literal.h"

#include "gundam/graph_type/large_graph.h"
#include "gundam/graph_type/large_graph2.h"
#include "gundam/graph_type/small_graph.h"
#include "gundam/graph_type/graph.h"

#include "gundam/graph_statistics/graph_basic_statistics.h"

#include "rule_discover/gar_discover/restriction.h"
#include "rule_discover/gar_discover/expand_tree.h"
#include "rule_discover/gar_discover/generate_tree.h"
#include "rule_discover/gar_discover/graph_package.h"

#include "util/log.h"

namespace grape {

using namespace _gar_discover;

/**
 * @brief Context for the parallel version of GarDiscover.
 *
 * @tparam FRAG_T
 */
template <typename FRAG_T>
class IncGarDiscoverContext : public VertexDataContext<FRAG_T, int64_t> {
 public:
  using oid_t = typename FRAG_T::oid_t;
  using vid_t = typename FRAG_T::vid_t;

  explicit IncGarDiscoverContext(const FRAG_T& fragment)
      : VertexDataContext<FRAG_T, int64_t>(fragment, true) {
    return;
  }

  void Init(ParallelMessageManager& messages, std::string yaml_file,
            int frag_num) {
    this->yaml_file_ = yaml_file;
    this->fid_ = this->fragment().fid();
    this->frag_num_ = frag_num;
#ifdef PROFILING
    preprocess_time = 0;
    exec_time = 0;
    postprocess_time = 0;
#endif
  }

  void Output(std::ostream& os) {
#ifdef PROFILING
    VLOG(2) << "preprocess_time: " << preprocess_time << "s.";
    VLOG(2) << "exec_time: " << exec_time << "s.";
    VLOG(2) << "postprocess_time: " << postprocess_time << "s.";
#endif
  }

#ifdef PROFILING
  double preprocess_time = 0;
  double exec_time = 0;
  double postprocess_time = 0;
#endif

  std::string yaml_file_;
  int fid_;
  int frag_num_;

  double time_limit_;
  double time_limit_per_supp_;

  // process num for each worker
  std::map<int, int> process_num_;

  using VertexIDType = uint64_t;
  using VertexLabelType = uint32_t;
  using VertexAttributeKeyType = std::string;

  using EdgeIDType = uint64_t;
  using EdgeLabelType = uint32_t;
  using EdgeAttributeKeyType = std::string;

  // using GraphPatternType =
  //     GUNDAM::LargeGraph2<VertexIDType, VertexLabelType, VertexAttributeKeyType,
  //                           EdgeIDType,   EdgeLabelType,   EdgeAttributeKeyType>;

  using GraphPatternType = GUNDAM::SmallGraph<VertexIDType, VertexLabelType, 
                                                EdgeIDType,   EdgeLabelType>;

  using DataGraphType =
      GUNDAM::LargeGraph<VertexIDType, VertexLabelType, VertexAttributeKeyType,
                           EdgeIDType,   EdgeLabelType,   EdgeAttributeKeyType>;
  using DataGraphLiteralType 
          = gar::LiteralInfo<GraphPatternType,
                                DataGraphType>;

  // using DataGraphType = GUNDAM::Graph<
  //     GUNDAM::SetVertexIDType<VertexIDType>,
  //     GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
  //     // GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
  //     GUNDAM::SetVertexLabelType<VertexLabelType>,
  //     GUNDAM::SetVertexLabelContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetVertexIDContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetVertexPtrContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetEdgeLabelContainerType<GUNDAM::ContainerType::Vector>,
  //     GUNDAM::SetVertexAttributeKeyType<VertexAttributeKeyType>,
  //     GUNDAM::SetEdgeIDType<EdgeIDType>,
  //     GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
  //     // GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
  //     GUNDAM::SetEdgeLabelType<EdgeLabelType>,
  //     GUNDAM::SetEdgeAttributeKeyType<EdgeAttributeKeyType>>;

  using GraphBasicStatisticsType = GUNDAM::GraphBasicStatistics<DataGraphType>;

  using RestrictionType = Restriction<GraphPatternType,
                                         DataGraphType>;

  using DataGraphVertexLabelType = typename GUNDAM::VertexLabel<DataGraphType>::type;
  using   DataGraphEdgeLabelType = typename GUNDAM::  EdgeLabel<DataGraphType>::type;

  using IncDataGraphPackage = IncGraphPackage<DataGraphType>;

  using RuleType = gar::GraphAssociationRule<GraphPatternType,
                                                DataGraphType>;
  using PatternVertexCounterType = typename GraphPatternType::VertexCounterType;
  using DataGraphVertexCounterType = typename DataGraphType::VertexCounterType;
  using GenerateTreeNodeType 
        = GenerateTreeNode<GraphPatternType,
                              DataGraphType>;
  using LiteralInfoType = gar::LiteralInfo<GraphPatternType,
                                              DataGraphType>;
  using GenerateTreeLevelType
      = GenerateTreeLevel<GraphPatternType,
                              DataGraphType>;

  bool RebuildLocalGARGenerateTree(std::unordered_map<int, int>& rule_id2parent_id,
                                std::unordered_map<int, DataGraphVertexCounterType>&
                                    rule_id2supp,
                                std::unordered_map<int, std::string>& rule_id2pivot_file,
                                std::unordered_map<int, double>& rule_id2probability,
                                int worker_id, int worker_num) {
    // input_gar_set must be ordered
    std::map<PatternVertexCounterType, std::map<unsigned, std::vector<unsigned>>> pattern_id2gar_set_vec;
    std::map<unsigned, unsigned> gar_id2pattern_id;

    for (unsigned gar_idx = 0; gar_idx < this->input_gar_set_.size(); gar_idx++) {
      PatternVertexCounterType vnum = this->input_gar_set_[gar_idx].pattern().CountVertex();
      if (gar_idx == 0) {
        unsigned pattern_id = gar_idx;
        pattern_id2gar_set_vec[vnum][pattern_id].push_back(gar_idx);
        gar_id2pattern_id[gar_idx] = pattern_id;;
        continue;
      }
      if (GUNDAM::SamePattern(this->input_gar_set_[gar_idx].pattern(), this->input_gar_set_[gar_idx - 1].pattern())) {
        unsigned pattern_id = gar_id2pattern_id[gar_idx - 1];
        pattern_id2gar_set_vec[vnum][pattern_id].push_back(gar_idx);
        gar_id2pattern_id[gar_idx] = pattern_id;
      } else {
        unsigned pattern_id = gar_idx;
        pattern_id2gar_set_vec[vnum][pattern_id].push_back(gar_idx);
        gar_id2pattern_id[gar_idx] = pattern_id;
      }
    }

    for (auto &[vnum, pattern_set] : pattern_id2gar_set_vec) {
      std::cout << "vnum " << vnum << " pattern num " << pattern_set.size() << std::endl;
    }


    std::map<unsigned, int> pattern2worker_id;
    for (auto &[level, pattern_id2gar_set] : pattern_id2gar_set_vec) {
        std::priority_queue<std::pair<unsigned, int>> gar_num_worker;
      for (int i = 0;i < worker_num; i++) {
        gar_num_worker.push(std::pair(0, i));
      }
      for (auto &[pattern_id, gar_set] : pattern_id2gar_set) {
        auto [least_worker_gar_num, least_worker_id] = gar_num_worker.top();
        gar_num_worker.pop();
        pattern2worker_id[pattern_id] = least_worker_id;
        //std::cout << "pattern id " <<  pattern_id << " least_worker_id " << least_worker_id << std::endl;
        gar_num_worker.push(std::make_pair(least_worker_gar_num
                                              + gar_set.size(), least_worker_id));
      }
    }
    
    std::map<unsigned, LiteralTreePtr> gar_id2literal_ptr;


    DataGraphVertexCounterType local_gar_num = 0;
    for (auto [pattern_id, current_worker_id] : pattern2worker_id) {
      if (current_worker_id != worker_id) {
        continue;
      }
      PatternVertexCounterType v_num = this->input_gar_set_[pattern_id].pattern().CountVertex();
      auto &current_pattern = this->input_gar_set_[pattern_id].pattern();
      GenerateTreeNodeType current_generate_node(pattern_id, current_pattern);

      typename GUNDAM::VertexID<LiteralTreeType>::type vertex_id_allocator = 1;
      typename GUNDAM::  EdgeID<LiteralTreeType>::type   edge_id_allocator = 0;

      for (unsigned gar_idx = pattern_id;
                    gar_idx < pattern_id + pattern_id2gar_set_vec[v_num][pattern_id].size();
                    gar_idx++) {
        if (rule_id2parent_id[gar_idx] == -1) {
          //std::cout << "y_literal root " << gar_idx << std::endl;
          auto &y_literal_set = this->input_gar_set_[gar_idx].y_literal_set();
          LiteralInfoType y_literal_info;

          for (auto &y_literal_ptr : y_literal_set) {
            y_literal_info = y_literal_ptr->info();
          }
          current_generate_node.AddLiteralTreeWithSupp(y_literal_info,
                                                      -1.0,
                                                      rule_id2supp[gar_idx],
                                                      gar_idx,
                                                      rule_id2pivot_file[gar_idx],
                                                      rule_id2probability[gar_idx]);
          auto& current_literal_tree = *(current_generate_node.LiteralTreesBackIter());
          LiteralTreePtr root_literal_handle = current_literal_tree.VertexBegin();
          local_gar_num++;
          gar_id2literal_ptr[gar_idx] = root_literal_handle;
          vertex_id_allocator = 1;
          edge_id_allocator = 0;
          continue;
        }

        auto parent_gar_idx = rule_id2parent_id[gar_idx];
        if (gar_id2literal_ptr.find(parent_gar_idx) == gar_id2literal_ptr.end()) {
          std::cout << "cant find a parent literal handle" << std::endl;
        }
        auto &parent_literal_handle = gar_id2literal_ptr[parent_gar_idx];

        auto &current_x_literal_set = this->input_gar_set_[gar_idx].x_literal_set();
        auto &parent_x_literal_set = this->input_gar_set_[parent_gar_idx].x_literal_set();

        bool found_same_literal = false;
        for (auto &current_x_literal_ptr : current_x_literal_set) {
          found_same_literal = false;
          for (auto &parent_x_literal_ptr : parent_x_literal_set) {
            if (current_x_literal_ptr->info() == parent_x_literal_ptr->info()) {
              found_same_literal = true;
              break;
            }
          }

          if (found_same_literal) {
            continue;
          }

          local_gar_num++;
          auto literal_tree_ptr = current_generate_node.LiteralTreesBackIter();
          auto& literal_tree = *literal_tree_ptr;
          auto [add_vertex_handle,
                add_vertex] = literal_tree.AddVertex(vertex_id_allocator++,
                                      parent_literal_handle->label() + 1);
          if (!add_vertex_handle) {
            std::cout << "fail to add vertex" << std::endl;
          }


          auto [ add_edge_handle,
                add_edge_ret ] = literal_tree.AddEdge(parent_literal_handle->id(), 
                                                            add_vertex_handle->id(),
                                                      kLiteralTreeDefaultEdgeLabel, 
                                                      edge_id_allocator++);
          if (!add_edge_ret) {
            std::cout << "fail to add edge" << std::endl;
          }

          // added successfully
          assert(add_edge_handle);
          assert(add_edge_ret);

          auto [ attr_handle,
                  attr_ret ] = add_vertex_handle->AddAttribute(kLiteralKey, 
                                                                current_x_literal_ptr->info());
                                                            
          // added successfully
          assert(attr_ret);

          auto [ supp_attr_handle,
                  supp_attr_ret ] = add_vertex_handle->AddAttribute(kSuppKey,
                                                                    rule_id2supp[gar_idx]);
                                                            
          assert(supp_attr_handle);
          assert(supp_attr_ret);
          //added successfully

          auto [ conf_attr_handle,
                conf_attr_ret ] = add_vertex_handle->AddAttribute(kConfKey, -1.0);
                
          // added successfully
          assert(conf_attr_ret);

          auto [pivot_file_handle,
                pivot_file_ret] = add_vertex_handle->AddAttribute(kMatchFileName, rule_id2pivot_file[gar_idx]);
          assert(pivot_file_handle);
          assert(pivot_file_ret);


          auto [origin_id_handle,
                origin_id_ret] = add_vertex_handle->AddAttribute(kOriginalGarID, gar_idx);
          assert(origin_id_handle);
          assert(origin_id_ret);

          auto [probability_handle,
                probability_ret] = add_vertex_handle->AddAttribute(kProbability, rule_id2probability[gar_idx]);
          assert(probability_handle);
          assert(probability_ret);

          gar_id2literal_ptr[gar_idx] = add_vertex_handle;

          break;
        }
        if (found_same_literal) {
          std::cout << "cant find the different literal " << gar_idx << std::endl;
        }
      }
      this->generate_tree_level_.AddGenerateTreeNode(std::move(current_generate_node));
    }
    std::cout << "worker id : " << worker_id << " local gar num " << local_gar_num << std::endl;
    return true;
  }



  RuleType PrepareGARandLiteralsToAdd(int gar_id, std::vector<DataGraphLiteralType>& eligible_literals) {
    RuleType gar = input_gar_set_[gar_id];

    int parent_id = rule_id2parent_id_[gar_id];
    if (parent_id == -1) {
      std::vector<DataGraphLiteralType> common_literal;
      for (auto &literal : gar.x_literal_set()) {
        common_literal.emplace_back(literal->info());
      }
      const auto& gar_pattern = input_gar_set_[gar_id].pattern();
      if (gar_id != 0) {
        for (int previous_gar_idx = gar_id - 1; previous_gar_idx >= 0; previous_gar_idx--) {
          const auto &previous_gar_pattern = input_gar_set_[previous_gar_idx].pattern();
          if (GUNDAM::SamePattern(gar_pattern, previous_gar_pattern)) {
            auto &same_pattern_gar = input_gar_set_[previous_gar_idx];
            for (auto &literal : same_pattern_gar.x_literal_set()) {
              bool found_same = false;
              for (unsigned literal_idx = 0; literal_idx < common_literal.size(); literal_idx++) {
                if (common_literal[literal_idx] == literal->info()) {
                  found_same = true;
                  break;
                }
              }
              if (!found_same) {
                common_literal.emplace_back(literal->info());
                eligible_literals.emplace_back(literal->info());
              }
            }
          }
        }
      }
      if (gar_id != input_gar_set_.size() - 1) {
        for (int after_gar_idx = gar_id - 1; after_gar_idx >= 0; after_gar_idx++) {
          const auto& after_gar_pattern = input_gar_set_[after_gar_idx].pattern();
          if (GUNDAM::SamePattern(gar_pattern, after_gar_pattern)) {
            auto &same_pattern_gar = input_gar_set_[after_gar_idx];
            for (auto &literal : same_pattern_gar.x_literal_set()) {
              bool found_same = false;
              for (unsigned literal_idx = 0; literal_idx < common_literal.size(); literal_idx++) {
                if (common_literal[literal_idx] == literal->info()) {
                  found_same = true;
                  break;
                }
              }
              if (!found_same) {
                common_literal.emplace_back(literal->info());
                eligible_literals.emplace_back(literal->info());
              }
            }
          }
        }
      }
      return gar;
    } else {
      std::vector<DataGraphLiteralType> common_literal;
      for (auto &literal : gar.x_literal_set()) {
        common_literal.emplace_back(literal->info());
      }
      for (auto &p : rule_id2parent_id_) {
        if ((p.second == parent_id) && (p.first != gar_id)) {
          auto &same_level_gar = input_gar_set_[p.first];
          for (auto &literal : same_level_gar.x_literal_set()) {
            bool found_same = false;
            for (unsigned literal_idx = 0; literal_idx < common_literal.size(); literal_idx++) {
              if (common_literal[literal_idx] == literal->info()) {
                found_same = true;
                break;
              }
            }
            if (!found_same) {
              common_literal.emplace_back(literal->info());
              eligible_literals.emplace_back(literal->info());
            }
          }
        }
      }
    }

    return gar;
  }


inline bool EligibleForExpandPattern(int gar_id) {
  bool eligible = true;
  if (rule_id2parent_id_[eligible] != -1) {
    return false;
  }
  const auto& gar_pattern = input_gar_set_[gar_id].pattern();
  if (gar_id != 0) {
    for (int pre_gar_id = gar_id - 1; pre_gar_id >= 0; pre_gar_id--) {
      const auto& pre_gar_pattern = input_gar_set_[pre_gar_id].pattern();
      if (!GUNDAM::SamePattern(gar_pattern, pre_gar_pattern)) {
        break;
      }
      if (rule_id2supp_[pre_gar_id] >= support_bound_) {
        eligible = false;
        break;
      }
    }
  }

  if (gar_id != input_gar_set_.size() - 1) {
    for (int after_gar_id = gar_id + 1; after_gar_id < input_gar_set_.size(); after_gar_id++) {
      const auto& after_gar_pattern = input_gar_set_[after_gar_id].pattern();
      if (!GUNDAM::SamePattern(gar_pattern, after_gar_pattern)) {
        break;
      }
      if (rule_id2supp_[after_gar_id] >= support_bound_) {
        eligible = false;
        break;
      }
    }
  }

  return eligible;
}


inline void CollectLiterals(int gar_id, std::vector<DataGraphLiteralType> &lhs_literals,
                                        std::vector<DataGraphLiteralType> &rhs_literals) {
  bool eligible = true;
  const auto &gar_pattern = input_gar_set_[gar_id].pattern();

  for (int target_gar_id = 0 ; target_gar_id < input_gar_set_.size(); target_gar_id++) {
    const auto& target_gar_pattern = input_gar_set_[target_gar_id].pattern();
    if (!GUNDAM::SamePattern(gar_pattern, target_gar_pattern)) {
      continue;
    }

    for (auto &literal : input_gar_set_[target_gar_id].y_literal_set()) {
      const auto &info = literal->info();
      bool has_same = false;
      for (unsigned i = 0; i < rhs_literals.size(); i++) {
        if (rhs_literals[i] == info) {
          has_same = true;
          break;
        }
      }
      if (!has_same) {
        rhs_literals.emplace_back(info);
      }
    }

    for (auto &literal : input_gar_set_[target_gar_id].x_literal_set()) {
      const auto &info = literal->info();
      bool has_same = false;
      for (unsigned i = 0; i < lhs_literals.size(); i++) {
        if (lhs_literals[i] == info) {
          has_same = true;
          break;
        }
      }
      if (!has_same) {
        lhs_literals.emplace_back(info);
      }
    }
  }
  return;
}



  std::vector<IncDataGraphPackage> data_graph_; // should have only one

  std::vector<RuleType> input_gar_set_;
  std::unordered_map<int, typename DataGraphType::VertexCounterType> rule_id2supp_;

  std::string knowledge_graph_v_file_,
              knowledge_graph_e_file_,
              er_file_;
  std::chrono::high_resolution_clock::time_point begin_time_;

  // GAR discover parameter
  int expand_round_;
  int j_; // literal tree depth

  // support bound
  uint64_t support_bound_;

  // to mark the current situation
  int current_round_;

  std::string output_gar_dir_;

  // GFD discover variable
  int root_pattern_max_edge_id_;

  bool with_index_, support_estimation_, match_gar_using_set_;

  double threshold_;

  // restrictions 
  RestrictionType restriction_;

  // hold the current level of the expand tree
  ExpandTreeLevel<GraphPatternType,
                     DataGraphType> expand_level_;

  std::ofstream time_log_file_;
  std::ofstream      log_file_;

  // to mark all expand node that needs to be processed
  std::vector<size_t> process_expand_node_id_set_;

  std::vector<std::vector<uint64_t>> pattern_match_count_;
  std::unordered_map<int, int> rule_id2parent_id_;
  // only be used in kExpandPatternFragID
  std::vector<std::pair<size_t, bool>> process_pattern_id_set_;

  GenerateTreeLevelType generate_tree_level_;

  // DataGraph Information
  GraphBasicStatisticsType origin_graph_basic_statistics_,
                          updated_graph_basic_statistics_;

  std::set<DataGraphEdgeLabelType> ml_edge_label_set_;

  std::set<std::tuple<DataGraphVertexLabelType,
                        DataGraphEdgeLabelType,
                      DataGraphVertexLabelType>> ml_edge_type_set_;
                  
  std::map<DataGraphEdgeLabelType,
           DataGraphEdgeLabelType> ml_to_normal_edge_label_,
                                   normal_to_ml_edge_label_;
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_INC_GAR_DISCOVER_GAR_DISCOVER_CONTEXT_H_
