#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GAR_DISCOVER_CONTEXT_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GAR_DISCOVER_CONTEXT_H_

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
class GarDiscoverContext : public VertexDataContext<FRAG_T, int64_t> {
 public:
  using oid_t = typename FRAG_T::oid_t;
  using vid_t = typename FRAG_T::vid_t;

  explicit GarDiscoverContext(const FRAG_T& fragment)
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

  // using DataGraphType =
  //     GUNDAM::LargeGraph2<VertexIDType, VertexLabelType, VertexAttributeKeyType,
  //                           EdgeIDType,   EdgeLabelType,   EdgeAttributeKeyType>;

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

  // using DataGraphType = GUNDAM::Graph<
  //     GUNDAM::SetVertexIDType<VertexIDType>,
  //     // GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
  //     GUNDAM::SetVertexAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
  //     GUNDAM::SetVertexLabelType<VertexLabelType>,
  //     GUNDAM::SetVertexLabelContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetVertexIDContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetVertexPtrContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetEdgeLabelContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetDecomposedEdgeContainerType<GUNDAM::ContainerType::Map>,
  //     GUNDAM::SetVertexAttributeKeyType<VertexAttributeKeyType>,
  //     GUNDAM::SetEdgeIDType<EdgeIDType>,
  //     // GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kGrouped>,
  //     GUNDAM::SetEdgeAttributeStoreType<GUNDAM::AttributeType::kSeparated>,
  //     GUNDAM::SetEdgeLabelType<EdgeLabelType>,
  //     GUNDAM::SetEdgeAttributeKeyType<EdgeAttributeKeyType>>;

  using GraphBasicStatisticsType = GUNDAM::GraphBasicStatistics<DataGraphType>;

  using RestrictionType = Restriction<GraphPatternType,
                                         DataGraphType>;

  using DataGraphVertexLabelType = typename GUNDAM::VertexLabel<DataGraphType>::type;
  using   DataGraphEdgeLabelType = typename GUNDAM::  EdgeLabel<DataGraphType>::type;

  using DataGraphPackage = GraphPackage<DataGraphType>;

  void ExportGraphSetStatistics(std::ostream& os) {
    os << "########################" << std::endl;
    for (size_t data_graph_idx = 0;
                data_graph_idx < this->data_graph_set_.size();
                data_graph_idx++) {
      os << "##  " << data_graph_idx 
         << "'th data graph: " << std::endl;

      GraphBasicStatisticsType graph_statistic(this->data_graph_set_[
                                                     data_graph_idx]
                                                    .data_graph());
      std::string graph_statistic_str = graph_statistic.ToString();
      os << graph_statistic_str;
    }
    return;
  }

  std::vector<DataGraphPackage> data_graph_set_;

  // GAR discover parameter
  int expand_round_;
  int j_; // literal tree depth

  // support bound
  uint64_t support_bound_;

  // to mark the current situation
  int current_round_;

  std::string output_gar_dir_,
              output_match_dir_;
  std::string knowledge_graph_v_file_,
              knowledge_graph_e_file_,
              er_file_;

  // GFD discover variable
  int root_pattern_max_edge_id_;

  // match gars with same pattern together
  bool match_gar_using_set_;
  bool store_match_;

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

  // only be used in kExpandPatternFragID
  std::vector<std::pair<size_t, bool>> process_pattern_id_set_;

  // DataGraph Information
  GraphBasicStatisticsType graph_basic_statistics_;

  std::set<DataGraphEdgeLabelType> ml_edge_label_set_;

  std::set<std::tuple<DataGraphVertexLabelType,
                        DataGraphEdgeLabelType,
                      DataGraphVertexLabelType>> ml_edge_type_set_;
                  
  std::map<DataGraphEdgeLabelType,
           DataGraphEdgeLabelType> ml_to_normal_edge_label_,
                                   normal_to_ml_edge_label_;
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GAR_DISCOVER_CONTEXT_H_
