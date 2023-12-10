#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GAR_DISCOVER_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GAR_DISCOVER_H_

#include <grape/grape.h>
#include <omp.h>

#include "../fragment2gundam.h"
#include "../fragment_graph.h"
#include "../fragment_graph_with_index.h"
#include "../timer.h"

#include "rule_discover/gar_discover_context.h"
#include "rule_discover/gar_discover/expand_literal.h"
#include "rule_discover/gar_discover/expand_pattern.h"
#include "rule_discover/gar_discover/export_gar.h"
#include "rule_discover/gar_discover/has_match.h"
#include "rule_discover/gar_discover/merge_pattern.h"
#include "rule_discover/gar_discover/serialize.h"

#include "grape/serialization/in_archive.h"
#include "grape/serialization/out_archive.h"
#include "gundam/io/csvgraph.h"
#include "gundam/algorithm/dp_iso.h"
#include "gundam/graph_type/large_graph2.h"
#include "gundam/graph_type/small_graph.h"

#include "gundam/tool/isolate_vertex.h"

#include "gar/literal_stand_alone_info.h"
#include "gar/same_gar.h"

#include "util/yaml_config.h"
#include "util/log.h"
#include "util/file.h"

namespace grape {

using namespace _gar_discover;

template <typename FRAG_T>
class GarDiscover : public ParallelAppBase<FRAG_T, GarDiscoverContext<FRAG_T>>,
                    public ParallelEngine {
 public:
  static constexpr LoadStrategy load_strategy = LoadStrategy::kLoadWholeGraph;

 private:
  using     VertexIDType = typename GarDiscoverContext<FRAG_T>::    VertexIDType;
  using  VertexLabelType = typename GarDiscoverContext<FRAG_T>:: VertexLabelType;
  using       EdgeIDType = typename GarDiscoverContext<FRAG_T>::      EdgeIDType;
  using    EdgeLabelType = typename GarDiscoverContext<FRAG_T>::   EdgeLabelType;
  using    DataGraphType = typename GarDiscoverContext<FRAG_T>::   DataGraphType;
  using GraphPatternType = typename GarDiscoverContext<FRAG_T>::GraphPatternType;
  
  using DataGraphLiteralType 
          = gar::LiteralInfo<GraphPatternType,
                                DataGraphType>;

  using VertexAttributeKeyType = typename GarDiscoverContext<FRAG_T>::VertexAttributeKeyType;

  using LiteralStandAloneInfoType 
 = gar::LiteralStandAloneInfo<GraphPatternType,
                                 DataGraphType>;

  using GarType = gar::GraphAssociationRule<GraphPatternType,
                                               DataGraphType>;

  static constexpr int kExpandPatternFragID = 0;
      
  // if pattern is not too many, preserve pattern with no legal
  // literal for pruning
  static const uint32_t kPatternTooMany = 5000;

  template <typename T>
  inline constexpr double CalTime(T& begin, T& end) {
    return double(std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                                        begin)
                      .count()) *
           std::chrono::microseconds::period::num /
           std::chrono::microseconds::period::den;
  }

  // from ordinary workers to kExpandPatternFragID
  //    info the kExpandPatternFragID how many proces the have
  static constexpr std::string_view kInfoProcessPrefix = "#info_process";
  // from kExpandPatternFragID to ordinary worker
  //    send pattern, new vertexes and possible literal
  static constexpr std::string_view kDeliverExpandNodesPrefix = "#deliver_expand_nodes";
  // from ordinary worker to kExpandPatternFragID
  //    update possible literals
  static constexpr std::string_view kUpdateExpandNodesPrefix = "#update_expand_nodes";

  // from kExpandPatternFragID to ordinary worker
  //    send pattern
  static constexpr std::string_view kDeliverPatternPrefix = "#deliver_pattern";
  // from ordinary worker to kExpandPatternFragID
  //    update pattern match count
  static constexpr std::string_view kUpdatePatternPrefix = "#update_pattern";

 public:
  // specialize the templated worker.
  INSTALL_PARALLEL_WORKER(GarDiscover<FRAG_T>, 
                          GarDiscoverContext<FRAG_T>,
                          FRAG_T)
  using vertex_t = typename fragment_t::vertex_t;

  /**
   * @brief Partial evaluation for GarDiscover.
   *
   * @param frag
   * @param ctx
   * @param messages
   */
  void PEval(const fragment_t& frag, 
                    context_t& ctx,
            message_manager_t& messages) {
    messages.InitChannels(thread_num());
    auto begin = timer();

    timer_next("load parameters");

    using namespace gar;

    std::ifstream config_file(ctx.yaml_file_);

    if (!config_file.is_open()) {
      util::Error("cannot open config file");
      return;
    }
    YAML::Node config = YAML::LoadFile(ctx.yaml_file_);
    
    // ##########################################
    // ##   load path(s) of data graph (set)   ##
    // ##########################################
    if (!config["DataGraphPath"]) {
      util::Error("cannot get data graph path");
      return;
    }
    YAML::Node data_graph_path_config = config["DataGraphPath"];

    std::vector<std::tuple<std::string,
                           std::string,
                           std::string>> data_graph_path_set;

    if (data_graph_path_config.IsSequence()) {
      data_graph_path_set.reserve(data_graph_path_config.size());
      for (int i = 0; i < data_graph_path_config.size(); i++) {
        auto [ data_graph_string, 
               data_graph_ret ] = util::GetDataGraphInfoFromYaml(data_graph_path_config[i]);
        if (!data_graph_ret) {
          util::Error("load data graph error!");
          return;
        }
        data_graph_path_set.emplace_back(data_graph_string);
      }
    }
    else {
      auto [ data_graph_string, 
             data_graph_ret ] = util::GetDataGraphInfoFromYaml(data_graph_path_config);
      if (!data_graph_ret) {
        util::Error("load data graph error!");
        return;
      }
      data_graph_path_set.emplace_back(data_graph_string);
    }
    
    // ##########################################
    // ##   discover algorithm configuration   ##
    // ##########################################
    if (!config["ExpandRound"]) {
      util::Error("cannot get expand round");
      return;
    }
    ctx. expand_round_ = config["ExpandRound"].as<int>();
    ctx.current_round_ = 0;
    
    // number of literals, contains Y
    if (!config["J"]) {
      util::Error("cannot get j");
      return;
    }
    ctx.j_ = config["J"].as<int>();

    // ####################################
    // ##   load specify edge type set   ##
    // ####################################
    using EdgeTypeType = std::tuple<VertexLabelType,
                                      EdgeLabelType,
                                    VertexLabelType>;

    std::set<EdgeTypeType> specified_edge_type_set;

    if (config["SpecifiedEdgeTypeSet"]) {
      YAML::Node specified_edge_type_set_config 
       = config["SpecifiedEdgeTypeSet"];
      for(int i=0; i < specified_edge_type_set_config.size(); ++i){
        YAML::Node specified_edge_type_config
                 = specified_edge_type_set_config[i];

        auto [edge_type, 
          get_edge_type_ret] = util::EdgeTypeFromYaml<DataGraphType>(specified_edge_type_config);

        if (!get_edge_type_ret) {
          // get edge type fail
          return;
        }

        auto [specified_edge_type_set_it,
              specified_edge_type_set_ret ]
            = specified_edge_type_set.emplace(edge_type);        
        if (! specified_edge_type_set_ret) {
          // add failed
          util::Error(std::string("specified duplicate edge type: <src_label, edge_label, dst_label>: ")
                    + "<" + std::to_string(std::get<0>(edge_type))
                    + "," + std::to_string(std::get<1>(edge_type))
                    + "," + std::to_string(std::get<2>(edge_type))
                    + ">");
          return;
        }
      }
      if (specified_edge_type_set.empty()) {
        util::Error("specified empty edge type");
        return;
      }
      util::Info("# Specified edge type number: " 
              + std::to_string(specified_edge_type_set.size())
              + " #");
    }

    // #####################################
    // ##   load specify edge label set   ##
    // #####################################
    std::set<EdgeLabelType> specified_edge_label_set;
    if (config["SpecifiedEdgeLabelSet"]) {
      YAML::Node specified_edge_label_set_config 
       = config["SpecifiedEdgeLabelSet"];
      for(int i=0; i < specified_edge_label_set_config.size(); ++i){
        EdgeLabelType specified_edge_label
                    = specified_edge_label_set_config[i].as<EdgeLabelType>();
        auto [specified_edge_label_set_it,
              specified_edge_label_set_ret ]
            = specified_edge_label_set.emplace(specified_edge_label);        
        if (! specified_edge_label_set_ret){
          util::Error("duplicate edge label: \t" 
                     + std::to_string(specified_edge_label));
          return;
        }
      }
      if (specified_edge_label_set.empty()){
        util::Error("specified empty edge label");
        return;
      }
    }

    if (config["ConstantFreqBound"]) {
      ctx.restriction_.set_constant_freq_bound(config["ConstantFreqBound"].as<float>());
      if (ctx.restriction_.constant_freq_bound() < 0.0
       || ctx.restriction_.constant_freq_bound() > 1.0){
        util::Error("illegal ConstantFreqBound: " + std::to_string(ctx.restriction_.constant_freq_bound()));
        return;
      }
    }

    if (config["KnowledgeGraphVFile"]) {
      ctx.knowledge_graph_v_file_ = config["KnowledgeGraphVFile"].as<std::string>();
    } else {
      ctx.knowledge_graph_v_file_ = "";
    }

    if (config["KnowledgeGraphEFile"]) {
      ctx.knowledge_graph_e_file_ = config["KnowledgeGraphEFile"].as<std::string>();
    } else {
      ctx.knowledge_graph_e_file_ = "";
    }

    if (config["ERFile"]) {
      ctx.er_file_ = config["ERFile"].as<std::string>();
    } else {
      ctx.er_file_ = "";
    }

    // ##########################################################
    // ##   load data graph (set) through the loaded path(s)   ##
    // ##########################################################
    for (const auto& [data_graph_v_file,
                      data_graph_e_file,
                   ml_literal_edge_file] : data_graph_path_set) {
      ctx.data_graph_set_.emplace_back(specified_edge_type_set,
                                       specified_edge_label_set,
                                       ctx.restriction_.constant_freq_bound(),
                                              data_graph_v_file,
                                              data_graph_e_file,
                                       ctx.graph_basic_statistics_,
                                       ctx.knowledge_graph_v_file_,
                                       ctx.knowledge_graph_e_file_,
                                       ctx.er_file_);
      if (!ctx.data_graph_set_.back().load_data_graph_success()) {
        util::Error("load data graph fail!");
        return;
      }
      ctx.data_graph_set_.back().GenerateGraphNec();
    }

    // #####################################################
    // ##  generate map from edge_label to ml_edge_label  ##
    // ##  as well as ml_edge_label to edge_label         ##
    // #####################################################
    EdgeLabelType max_edge_label = 0,
      minimal_negtive_edge_label = 0;
    for (const auto& [edge_label,
                      edge_label_counter]
                : ctx.graph_basic_statistics_
                     .edge_label_counter()) {
      if (edge_label < 0) {
        if (minimal_negtive_edge_label > edge_label) {
          minimal_negtive_edge_label = edge_label;
        } 
        continue;
      }
      if (max_edge_label < edge_label) {
        max_edge_label = edge_label;
      }
    }
    EdgeLabelType edge_label_offset = max_edge_label 
                        - minimal_negtive_edge_label;

    for (const auto& [edge_label, edge_label_counter] 
                            : ctx.graph_basic_statistics_
                                 .edge_label_counter()) {
      assert(ctx.ml_to_normal_edge_label_.find(edge_label + edge_label_offset)
          == ctx.ml_to_normal_edge_label_.end());
      assert(ctx.normal_to_ml_edge_label_.find(edge_label)
          == ctx.normal_to_ml_edge_label_.end());
      ctx.ml_to_normal_edge_label_.emplace(edge_label + edge_label_offset, edge_label);
      ctx.normal_to_ml_edge_label_.emplace(edge_label, edge_label + edge_label_offset);

      ctx.ml_edge_label_set_.emplace(edge_label + edge_label_offset);
    }

    for (const auto& [edge_type, edge_type_counter] 
                           : ctx.graph_basic_statistics_
                                .edge_type_counter()) {
      auto [ ml_edge_type_set_it,
             ml_edge_type_set_ret ]
       = ctx.ml_edge_type_set_.emplace(std::get<0>(edge_type),  // src_label
                                       std::get<1>(edge_type) + edge_label_offset,
                                       std::get<2>(edge_type)); // dst_label
      assert(ml_edge_type_set_ret);
    }
    
    // ##########################
    // ##   load ml edge set   ##
    // ##########################
    for (int data_graph_idx = 0; 
             data_graph_idx < data_graph_path_set.size(); 
             data_graph_idx++) {

      const std::string& ml_e_file = std::get<2>(data_graph_path_set[data_graph_idx]);
    
      if (ml_e_file == "") {
        // does not have ml edge for this data graph
        continue;
      }
      if (!ctx.data_graph_set_[data_graph_idx].AddMlEdge(
                ml_e_file,
                ctx.normal_to_ml_edge_label_)) {
        util::Error("add ml edge fail!");
        return;
      }
    }

    if (config["PatternVertexLimit"]) {
      ctx.restriction_.set_pattern_vertex_limit(config["PatternVertexLimit"].as<int>());
    }

    if (config["DiameterLimit"]) {
      ctx.restriction_.set_diameter_limit(config["DiameterLimit"].as<int>());
    }

    ctx.support_bound_ = 1;
    if (config["SupportBound"]) {
      ctx.support_bound_ = config["SupportBound"].as<decltype(ctx.support_bound_)>();
    }

    if (config["ConfidenceBound"]) {
      assert(!ctx.restriction_.specified_confidence_bound());
      ctx.restriction_.set_confidence_bound(config["ConfidenceBound"].as<float>());

      if (ctx.restriction_.confidence_bound() < 0.0
       || ctx.restriction_.confidence_bound() > 1.0) {
        util::Error("illegal ConfidenceBound: " + std::to_string(ctx.restriction_.confidence_bound()));
        return;
      }
      assert(ctx.restriction_.specified_confidence_bound());
    }

    if (!config["OutputGarDir"]) {
      util::Error("cannot get output gar dir");
      return;
    }
    ctx.output_gar_dir_ = config["OutputGarDir"].as<std::string>();

    util::Mkdir(ctx.output_gar_dir_);

    // ########################################
    // ##   load specified Rhs literal set   ##
    // ########################################
    if (config["SpecifiedRhsLiteralSet"]) {
      YAML::Node specified_rhs_literal_set_config 
       = config["SpecifiedRhsLiteralSet"];
      for (int specified_rhs_literal_set_config_idx = 0; 
               specified_rhs_literal_set_config_idx 
             < specified_rhs_literal_set_config.size(); 
               specified_rhs_literal_set_config_idx++) {
        const auto& specified_rhs_literal_config
                  = specified_rhs_literal_set_config[specified_rhs_literal_set_config_idx];
        auto [ stand_alone_literal_info,
               stand_alone_literal_info_ret ]
          = util::LiteralStandAloneInfoFromYaml<GraphPatternType,
                                                   DataGraphType>(specified_rhs_literal_config);
        if (!stand_alone_literal_info_ret) {
          util::Error("the " + std::to_string(specified_rhs_literal_set_config_idx) 
                             + "'th literal is illegal");
          return;
        }
        if (!specified_rhs_literal_config["SuppFile"]) {
          // does no specifed support file
          ctx.restriction_.AddSpecifiedRhsLiteral(stand_alone_literal_info);
          continue;
        }
        if (ctx.data_graph_set_.size() > 1) {
          util::Error("support set only allowed to be specified on one data graph");
          return;
        }
        // ################################
        // ##   specified support file   ##
        // ################################
        const size_t kLiteralVertexNum = stand_alone_literal_info.VertexNum();
        if (kLiteralVertexNum != 2) {
          util::Error("support only literals contains two vertex now!");
          return;
        }
        if (stand_alone_literal_info.literal_type() != gar::LiteralType::kEdgeLiteral) {
          // variable literal has problem in simulation with specified support set
          // in HasMatch
          util::Error("support only edge literals contains two vertex now!");
          return;
        }

        const std::string kLiteralSuppFile 
           = specified_rhs_literal_config["SuppFile"].as<std::string>();

        std::ifstream literal_supp_file(kLiteralSuppFile);

        if (!literal_supp_file.good()) {
          util::Error ( "literal_supp_file of the " 
                      + std::to_string(specified_rhs_literal_set_config_idx) 
                      + "'th literal: " + kLiteralSuppFile
                      + " is not good! ");
          return;
        }

        std::vector<
        std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>> literal_supp_set;

        while (literal_supp_file) {
          std::string s;
          if (!std::getline( literal_supp_file, s ))
            break;

          std::istringstream ss( s );
          std::vector <std::string> record;
          std::string buf; 
          while (ss >> buf)
            record.emplace_back(buf);

          std::vector <VertexIDType> support_id_set;
          for (const auto& id_str : record) {
            support_id_set.emplace_back(std::stoi(id_str));
          }
          if (kLiteralVertexNum != support_id_set.size()) {
            if ((kLiteralVertexNum + 1) != support_id_set.size()
             || (support_id_set[kLiteralVertexNum] != 0
              && support_id_set[kLiteralVertexNum] != 1)) {
              util::Error("support set miss match in line: " 
                         + std::to_string(literal_supp_set.size())
                         + " of the " 
                         + std::to_string(specified_rhs_literal_set_config_idx) 
                         + "'th literal");
              return;
            }
          }
          if (kLiteralVertexNum + 1 == support_id_set.size()) {
            assert(support_id_set[kLiteralVertexNum] == 0
                || support_id_set[kLiteralVertexNum] == 1);
            support_id_set.pop_back();
          }
          std::vector <typename GUNDAM::VertexHandle<DataGraphType>::type> support;
          support.reserve(support_id_set.size());
          // support set only allowed to be specified on one data graph
          assert(ctx.data_graph_set_.size() == 1);
          for (const auto& vertex_id : support_id_set) {
            auto vertex_handle = ctx.data_graph_set_.at(0)
                                    .data_graph().FindVertex(vertex_id);
            if (!vertex_handle) {
              util::Error("vertex: " + std::to_string(vertex_id) 
                        + " cannot be found in data graph!");
              return;
            }
            support.emplace_back(vertex_handle);
          }
          assert(support.size() == support_id_set.size());
          literal_supp_set.emplace_back(support);
        }
        if (literal_supp_set.empty()) {
          util::Error("Specified empty support set for the "
                     + std::to_string(specified_rhs_literal_set_config_idx) 
                     + "'th literal");
          return;
        }
        std::sort(literal_supp_set.begin(),
                  literal_supp_set.end());
        ctx.restriction_.AddSpecifiedRhsLiteral(stand_alone_literal_info,
                                                literal_supp_set);
      }
      util::Info("specified_rhs_literal_set size: "
                 + std::to_string(ctx.restriction_.specified_rhs_literal_set().size()));
    }

    ctx.time_limit_ = -1.0;
    if (config["TimeLimit"]) {
      ctx.time_limit_ = config["TimeLimit"].as<double>();
    }
    
    ctx.time_limit_per_supp_ = -1.0;
    if (config["TimeLimitPerSupp"]) {
      ctx.time_limit_per_supp_ = config["TimeLimitPerSupp"].as<double>();
    }

    if (config["TimeLogFile"]) {
      assert(!ctx.time_log_file_.is_open());
      const std::string kTimeLogFileName = config["TimeLogFile"].as<std::string>();
      ctx.time_log_file_.open(kTimeLogFileName, std::ios::app);
      if (!ctx.time_log_file_.is_open()) {
        util::Error("open time log file failed: " + kTimeLogFileName);
        return;
      }
      ctx.time_log_file_ << "time log of: " << ctx.yaml_file_ << std::endl;
    }

    if (config["LogFile"]) {
      /* ######################## *
       * ##  for detailed log  ## *
       * ######################## */
      assert(!ctx.log_file_.is_open());
      const std::string kLogFileName = config["LogFile"].as<std::string>();
      ctx.log_file_.open(kLogFileName, std::ios::app);
      if (!ctx.log_file_.is_open()) {
        util::Error("open log file failed: " + kLogFileName);
        return;
      }
      ctx.log_file_ << "log of: " << ctx.yaml_file_ << std::endl;
    }

    // ################################################
    // ##   load types of literal to be considered   ##
    // ################################################
    if (config["LiteralTypes"]) {
      ctx.restriction_.set_consider_no_literal();
      YAML::Node literal_types_config = config["LiteralTypes"];
      for (int i = 0; i < literal_types_config.size(); i++) {
        const std::string literal_type 
                        = literal_types_config[i].as<std::string>();
        if (!ctx.restriction_.ConsiderLiteral(literal_type)) {
          // illegal literal
          return;
        }
      }
    }

    // #######################################################
    // ##   load types of literal to be considered in LHS   ##
    // #######################################################
    if (config["LhsLiteralTypes"]) {
      ctx.restriction_.set_consider_no_literal_in_lhs();
      YAML::Node lhs_literal_types_config = config["LhsLiteralTypes"];
      for (int i = 0; i < lhs_literal_types_config.size(); i++) {
        const std::string lhs_literal_type 
                        = lhs_literal_types_config[i].as<std::string>();
        if (!ctx.restriction_.ConsiderLiteralInLhs(lhs_literal_type)) {
          // illegal literal
          return;
        }
      }
    }

    // ##########################################
    // ##   load graph pattern configuration   ##
    // ##########################################
    if (config["Pattern"]) {
      YAML::Node pattern_config = config["Pattern"];
      
      std::string pattern_type_str = pattern_config["Type"].as<std::string>();

      if (!ctx.restriction_.SetPatternType(pattern_type_str)) {
        util::Error("PatternType error: "+ pattern_type_str);
        return;
      }

      if (ctx.restriction_.pattern_is_star()) {
        YAML::Node star_config = pattern_config;
        assert(!ctx.restriction_.specified_path_num_limit());
        if (star_config["PathNumLimit"]) {
          const auto kPathNumLimit = star_config["PathNumLimit"].as<int>();
          if (kPathNumLimit < 1) {
            util::Error("illegal PathNumLimit for gcr: " 
               + std::to_string(kPathNumLimit));
            return;
          }
          ctx.restriction_.set_path_num_limit(kPathNumLimit);
          assert(ctx.restriction_.path_num_limit() == kPathNumLimit);
          assert(ctx.restriction_.specified_path_num_limit());
        }
        assert(!ctx.restriction_.specified_path_length_limit());
        if (star_config["PathLengthLimit"]) {
          const auto kPathLengthLimit = star_config["PathLengthLimit"].as<int>();
          if (kPathLengthLimit < 1) {
            util::Error("illegal PathLengthLimit for gcr: " 
               + std::to_string(kPathLengthLimit));
            return;
          }
          ctx.restriction_.set_path_length_limit(kPathLengthLimit);
          assert(ctx.restriction_.path_length_limit() == kPathLengthLimit);
          assert(ctx.restriction_.specified_path_length_limit());
        }
        if (star_config["BidirectionalPath"]) {
          const auto kBidirectionalPath = star_config["BidirectionalPath"].as<bool>();
          ctx.restriction_.set_bidirectional_path(kBidirectionalPath);
          assert(ctx.restriction_.bidirectional_path() == kBidirectionalPath);
        }
      }
      else if (ctx.restriction_.pattern_is_tree()) {
        util::Error("does not support configure for tree");
        return;
      }
      else if (ctx.restriction_.pattern_is_link()) {
        util::Error("does not support configure for link");
        return;
      }
    }
    
    // #####################################
    // ##   load rule types to discover   ##
    // #####################################
    if (config["Rule"]) {
      YAML::Node rule_config = config["Rule"];
      const std::string rule_type 
                      = rule_config["Type"].as<std::string>();
      if (!ctx.restriction_.SpecifyRuleType(rule_type)) {
        // specify rule type fail
        return;
      }
      if (ctx.restriction_.gcr()) {
        if (rule_config["Path"]) {
          YAML::Node path_config = rule_config["Path"];
          for (int i = 0; i < path_config.size(); i++) {
            auto [ path_string, 
                   path_ret ] = util::GetDataGraphInfoFromYaml(path_config[i]);
            if (!path_ret) {
              util::Error("load path error!");
              return;
            }
            const auto& v_file_path = std::get<0>(path_string);
            const auto& e_file_path = std::get<1>(path_string);
            GraphPatternType path;
            if (GUNDAM::ReadCSVGraph(path, v_file_path,
                                           e_file_path) < 0) {
              util::Error("load path fail!");
              return;
            }
            ctx.restriction_.AddSpecifiedPath(path);
          }
        }
        if (rule_config["Restrictions"]) {
          YAML::Node restrictions_config = rule_config["Restrictions"];
          for (int i = 0; i < restrictions_config.size(); i++) {
            const std::string restriction 
                            = restrictions_config[i].as<std::string>();
            if (!ctx.restriction_.AddRestriction(restriction)) {
              // add restriction failed!
              return;
            }
          }
        }
        if (ctx.restriction_.central_to_central()) {
          if (!ctx.restriction_.discretized_graph()) {
            util::Error("GCR: only support central_to_central in discretized_graph");
            return;
          }
          for (auto& data_graph : ctx.data_graph_set_) {
            data_graph.GenerateCentralSet();
          }
        }
        if (rule_config["CentralVertexSample"]) {
          if (!ctx.restriction_.central_to_central()) {
            util::Error("specified sample ratio but not specified central_to_central");
            return;
          }
          const double kSampleRatio = rule_config["CentralVertexSample"].as<double>();
          if (kSampleRatio <= 0 || kSampleRatio  >= 1) {
            util::Error("Illegal sample ratio: " + std::to_string(kSampleRatio));
            return;
          }
          for (auto& data_graph : ctx.data_graph_set_) {
            data_graph.SampleCentralSet(kSampleRatio);
          }
        }
      }
    }
    else {
      ctx.restriction_.SpecifyRuleType("gar");
    }
    assert(ctx.restriction_.gar()
        || ctx.restriction_.gcr()
        || ctx.restriction_.horn_rule());

    // #######################################
    // ##load whether matching gars as a set##
    // #######################################

    if (config["MatchGARUsingSet"]) {
      ctx.match_gar_using_set_ = config["MatchGARUsingSet"].as<bool>();
      if (ctx.match_gar_using_set_) {
        std::cout << "####match gar using set#####" << std::endl;
      } else {
        std::cout << "####match gar one by one####" << std::endl;
      }
    } else {
      std::cout << "Default : ####match gar one by one####" << std::endl;
      ctx.match_gar_using_set_ = false;
    }

    if (config["StoreMatch"]) {
      ctx.store_match_ = config["StoreMatch"].as<bool>();
      if (ctx.store_match_) {
        std::cout << "######gar discover will store match######" << std::endl;
      } else {
        std::cout << "######gar discover does not store match######" << std::endl;
      }
    } else {
      std::cout << "Default : ###gar discover does not store match####" << std::endl;
      ctx.store_match_ = false;
    }

    if (config["OutputMatchDir"]) {
      ctx.output_match_dir_ = config["OutputMatchDir"].as<std::string>();
      std::cout << "output match dir " << ctx.output_match_dir_ << std::endl;
      util::Mkdir(ctx.output_match_dir_);
    } else {
      if (ctx.store_match_) {
        util::Error("should store match but have no output match dir");
        return;
      }
    }

    // #################################
    // ##   load other restrictions   ##
    // #################################
    if (config["Restrictions"]) {
      YAML::Node restrictions_config = config["Restrictions"];
      for (int i = 0; i < restrictions_config.size(); i++) {
        const std::string restriction 
                        = restrictions_config[i].as<std::string>();
        if (!ctx.restriction_.AddRestriction(restriction)) {
          // add restriction failed!
          return;
        }
      }
    }

    if (ctx.fid_ == kExpandPatternFragID
     && ctx.log_file_.is_open()) {
      // export log
      ctx.ExportGraphSetStatistics(ctx.log_file_);
    }

    if (ctx.fid_ != kExpandPatternFragID) {
      // has load all parameters and has finished
      // all configuration
      std::string msg(kInfoProcessPrefix);
      msg += " " + std::to_string(ctx.fid_)
           + " " + std::to_string(omp_get_num_procs());
      auto& channel_0 = messages.Channels()[0];
      channel_0.SendToFragment(kExpandPatternFragID, msg);
      return;
    }

    // is PatternExpand node, patterns are expanded in this node
    // make other preparations
    std::vector<GraphPatternType> root_patterns;

    std::set<std::pair<VertexLabelType,
                       VertexLabelType>> considered_label_pair_set;
    for (const auto& [vertex_label, 
                      vertex_label_counter]
                : ctx.graph_basic_statistics_
                     .vertex_label_counter()) {
      GraphPatternType root_pattern;
      root_pattern.AddVertex(0, vertex_label);
      if (!ctx.restriction_.gcr()) {
        root_patterns.emplace_back(std::move(root_pattern));
        continue;
      }
      for (const auto& [other_vertex_label, 
                        other_vertex_label_counter]
                  : ctx.graph_basic_statistics_
                       .vertex_label_counter()) {
        auto [considered_label_pair_set_it,
              considered_label_pair_set_ret] 
            = considered_label_pair_set.emplace(std::min(vertex_label, other_vertex_label),
                                                std::max(vertex_label, other_vertex_label));
        if (!considered_label_pair_set_ret) {
          // has been considered
          continue;
        }
        GraphPatternType two_vertexes_root_pattern(root_pattern);
        two_vertexes_root_pattern.AddVertex(1, other_vertex_label);
        root_patterns.emplace_back(std::move(two_vertexes_root_pattern));
      }
    }

    timer_next("run gar_discover level 0");

    ctx.expand_level_.clear();

    // expand from the root graph pattern(s)
    ctx.root_pattern_max_edge_id_ = 0;
    int pattern_counter = 0;
    for (const auto& root_pattern : root_patterns) {
      assert(ctx.root_pattern_max_edge_id_ == 0 
          || ctx.root_pattern_max_edge_id_ == root_pattern.CountEdge());
      // does not support specified pattern yet
      assert(root_pattern.CountEdge() == 0);
      std::vector<VertexIDType> first_considered_ids;
      // all vertexes are considered at the first time
      for (auto vertex_cit = root_pattern.VertexBegin();
               !vertex_cit.IsDone(); 
                vertex_cit++) {
        first_considered_ids.emplace_back(vertex_cit->id());
      }
      auto& expand_tree_node 
          = ctx.expand_level_.AddExpandTreeNode(
                                pattern_counter++, 
                                root_pattern, 
                                first_considered_ids,
                                // possible rhs literals
                                std::vector<DataGraphLiteralType>(),
                                // all lhs literals
                                std::vector<DataGraphLiteralType>()
                              );
      // initilize the possible literal of the root pattern
      std::vector<DataGraphLiteralType> possible_literals;
      ExpandRhsLiteral(expand_tree_node, 
                       ctx.graph_basic_statistics_,
                       ctx.ml_edge_label_set_,
                       ctx.ml_edge_type_set_,
                       ctx.ml_to_normal_edge_label_,
                       ctx.restriction_,
                       possible_literals);

      std::string root_pattern_str;
      root_pattern_str << root_pattern;
      util::Info("root_pattern_str: "
                + root_pattern_str);
                       
      util::Info("possible_literals.size(): "
                + std::to_string(possible_literals.size()));
      for (const auto& literal : possible_literals) {
        std::string literal_str;
        literal_str << literal;
        util::Info(literal_str);
        // does not add to rhs literal
        if (literal.literal_type() != gar::LiteralType::kMlLiteral) {
          // ml literal won't be added into rhs
          expand_tree_node.AddRhsLiteral(literal);
        }
        if (ctx.j_ == 1) {
          continue;
        }
        if (literal.literal_type() != gar::LiteralType::kEdgeLiteral) {
          // edge literal won't be added into lhs
          if (ctx.restriction_.LiteralTypeConsideredInLhs(literal.literal_type())) {
            expand_tree_node.AddLhsLiteral(literal);
          }
        }
      }
      assert(expand_tree_node.CountNewVertexes() != 0);
    }

    assert(ctx.expand_level_.size() != 0);
    util::Info("root patterns size: " + std::to_string(ctx.expand_level_.size()));
      
    std::string msg(kInfoProcessPrefix);
    msg += " " + std::to_string(ctx.fid_)
         + " " + std::to_string(omp_get_num_procs());
    auto& channel_0 = messages.Channels()[0];
    channel_0.SendToFragment(kExpandPatternFragID, msg);

    // #ifdef PROFILING
    // ctx.exec_time -= GetCurrentTime();
    // #endif

    // #ifdef PROFILING
    // ctx.exec_time += GetCurrentTime();
    // ctx.postprocess_time -= GetCurrentTime();
    // #endif
    // // messages.ForceContinue();

    // #ifdef PROFILING
    // ctx.postprocess_time += GetCurrentTime();
    // #endif
  }

  /**
   * @brief Incremental evaluation for GarDiscover.
   *
   * @param frag
   * @param ctx
   * @param messages
   */
  void IncEval(const fragment_t& frag, 
                      context_t& ctx,
              message_manager_t& messages) {
    // auto inner_vertices = frag.InnerVertices();

    using ExpandTreeNodeType 
        = ExpandTreeNode<GraphPatternType,
                            DataGraphType>;

    using ExpandTreeLevelType
        = ExpandTreeLevel<GraphPatternType,
                             DataGraphType>;
                      
    using GenerateTreeNodeType 
        = GenerateTreeNode<GraphPatternType,
                              DataGraphType>;
 
    using GenerateTreeLevelType
        = GenerateTreeLevel<GraphPatternType,
                               DataGraphType>;

    auto& channels = messages.Channels();
    // using VertexIDType = typename fragment_t::vid_t;

    bool receive_info_process        = false,
         receive_deliver_pattern     = false,
         receive_update_pattern      = false,
         receive_deliver_expand_node = false,
         receive_update_expand_node  = false,
         receive_message             = false;

    std::vector<std::tuple<GraphPatternType, size_t, bool>> delivered_patterns;

    assert(ctx.fid_ == kExpandPatternFragID
        || ctx.expand_level_.empty());

    util::Info("omp_get_num_procs: "
              + std::to_string(omp_get_num_procs()));

    messages.ParallelProcess<std::string>(
        // thread_num(),
        1, [&ctx, 
            &delivered_patterns,
            &receive_info_process,
            &receive_deliver_pattern,
            &receive_update_pattern,
            &receive_deliver_expand_node, 
            &receive_update_expand_node, 
            &receive_message](int tid, std::string msg) {
          util::Info("fid: " + std::to_string(ctx.fid_) 
                  + " receive message: " + msg);

          auto res_update_expand_node 
             = std::mismatch(kUpdateExpandNodesPrefix.begin(),
                             kUpdateExpandNodesPrefix. end (), msg.begin());
          auto res_deliver_expand_node
             = std::mismatch(kDeliverExpandNodesPrefix.begin(),
                             kDeliverExpandNodesPrefix. end (), msg.begin());
          auto res_update_pattern 
             = std::mismatch(kUpdatePatternPrefix.begin(),
                             kUpdatePatternPrefix. end (), msg.begin());
          auto res_deliver_pattern =
              std::mismatch(kDeliverPatternPrefix.begin(),
                            kDeliverPatternPrefix. end (), msg.begin());
          auto res_info =     
              std::mismatch(kInfoProcessPrefix.begin(),
                            kInfoProcessPrefix. end (), msg.begin());
          receive_message = true;
          if (res_update_expand_node.first == kUpdateExpandNodesPrefix.end()) {
            msg = msg.substr(kUpdateExpandNodesPrefix.size());
            // kUpdateExpandNodesPrefix is the prefix of msg.
            receive_update_expand_node = true;
            assert(ctx.fid_ == kExpandPatternFragID);
            // ###############################################################
            // ##  receive the message from the ordinary workers after the  ##
            // ##  horizontal spawning, update the information in current   ##
            // ##  level expand tree                                        ##
            // ###############################################################
            while (!msg.empty()) {
              UpdateLiteral<DataGraphLiteralType> update_literal;
              msg >> update_literal;
              ctx.expand_level_.node(update_literal.id())
                               .rhs_literals()
                               .swap(update_literal.rhs_literals());
              ctx.expand_level_.node(update_literal.id())
                               .lhs_literals()
                               .swap(update_literal.lhs_literals());
            }
          } else if (res_deliver_expand_node.first == kDeliverExpandNodesPrefix.end()) {
            msg = msg.substr(kDeliverExpandNodesPrefix.size());
            // kDeliverExpandNodesPrefix is the prefix of msg.
            receive_deliver_expand_node = true;
            // #########################################################################
            // ##  receive the expand nodes delievered from the kExpandPatternFragID  ##
            // ##  for horizontal spawing to evaluate the literals on it              ##
            // #########################################################################
            while (!msg.empty()) {
              ExpandTreeNodeType temp_expand_tree_node;
              msg >> temp_expand_tree_node;
              ctx.expand_level_.AddExpandTreeNode(temp_expand_tree_node);
            }
            util::Info("#delivered expand node#");
            util::Info("ctx.expand_level_.size(): "
                      + std::to_string(ctx.expand_level_.size()));
            util::Info("expand node number: " 
                      + std::to_string(
                          (ctx.fid_ == kExpandPatternFragID) ?
                           ctx.process_expand_node_id_set_.size()
                         : ctx.expand_level_.size() ));
          } else if (res_deliver_pattern.first == kDeliverPatternPrefix.end()) {
            msg = msg.substr(kDeliverPatternPrefix.size());
            // kDeliverPatternPrefix is the prefix of msg.
            receive_deliver_pattern = true;
            // #####################################################################
            // ##  receive the patterns delievered from the kExpandPatternFragID  ##
            // ##  to evaluate the match count                                    ##
            // #####################################################################
            while (!msg.empty()) {
              std::tuple<GraphPatternType, size_t, bool> temp_delivered_pattern;
              msg >> temp_delivered_pattern;
              delivered_patterns.emplace_back(std::move(temp_delivered_pattern));
            }
          } else if (res_update_pattern.first == kUpdatePatternPrefix.end()) {
            msg = msg.substr(kUpdatePatternPrefix.size());
            // kDeliverPatternPrefix is the prefix of msg.
            receive_update_pattern = true;
            assert(ctx.fid_ == kExpandPatternFragID);
            // ###########################################################
            // ##  receive the message from the ordinary workers after  ##
            // ##  evaluating the match count                           ##
            // ###########################################################
            std::stringstream ss;
            if (msg != "") {
              std::istringstream ss( msg );
              while (ss) {
                std::string pattern_id_str,
                           match_count_str,
                              graph_id_str;
                ss >>  pattern_id_str;
                if (pattern_id_str == "") {
                  break;
                }
                ss >> match_count_str;
                if (match_count_str == "") {
                  util::Error("# Error! cannot get match_count_str! #");
                  break;
                }
                ss >> graph_id_str;
                if (graph_id_str == "") {
                  util::Error("# Error! cannot get graph_id_str! #");
                  break;
                }
                util::Info("\n  pattern_id_str:" +  pattern_id_str + "#\n"
                           + " match_count_str:" + match_count_str + "#\n"
                           + "    graph_id_str:" +    graph_id_str + "#\n");
                size_t pattern_id = std::stoi(pattern_id_str);
                if (pattern_id >= ctx.pattern_match_count_.size()
                 || pattern_id < 0) {
                  util::Error("pattern_id >= ctx.pattern_match_count_.size()!\n");
                  util::Error("pattern_id: " + std::to_string(pattern_id)
                           + " ctx.pattern_match_count_.size(): " + std::to_string(ctx.pattern_match_count_.size()));
                }
                size_t graph_id = std::stoi(graph_id_str);
                if (graph_id >= ctx.pattern_match_count_[pattern_id].size()
                 || graph_id < 0) {
                  util::Error("graph_id >= ctx.pattern_match_count_[pattern_id].size()!");
                  util::Error("graph_id: " + std::to_string(graph_id)
                           + " ctx.pattern_match_count_[pattern_id].size(): " + std::to_string(ctx.pattern_match_count_[pattern_id].size()));
                }
                uint64_t match_count = std::stoi(match_count_str);
                ctx.pattern_match_count_[pattern_id][graph_id] = match_count;
              }
            }
            assert(ctx.expand_level_.size() == ctx.pattern_match_count_.size());
          } else if (res_info.first == kInfoProcessPrefix.end()) {
            msg = msg.substr(kInfoProcessPrefix.size());
            receive_info_process = true;
            assert(ctx.fid_ == kExpandPatternFragID);
            // ###############################################################
            // ##  receive the number of processors of the ordinary workers  #
            // ###############################################################
            std::string fid_str, proc_num_str;
            std::stringstream ss;
            ss << msg;
            ss >> fid_str;
            ss >> proc_num_str;
            util::Info(" fid: " + fid_str + " proc_num: " + proc_num_str);
            ctx.process_num_.emplace(std::stoi(fid_str), 
                                     std::stoi(proc_num_str));
          } else {
            // unknown message type
            assert(false);
          }
        });

    if (!receive_message) {
      return;
    }

    if (ctx.fid_ == kExpandPatternFragID) {
      // only the kExpandPatternFragID can receive feedback
      // from other workers

      if (receive_update_expand_node 
       || receive_info_process) {
        // ##################################################
        // ##  receive the infomation from the workers to  ##
        // ##  update the expand node for the next round   ##
        // ##################################################

        if (ctx.current_round_ > 0) {
          if (ctx.time_log_file_.is_open()) {
            ctx.time_log_file_ << timer_now() << std::endl;
          }
          timer_next("run gar_discover level " + std::to_string(ctx.current_round_));
        }
        
        util::Output( " current round: " + std::to_string(ctx.current_round_)
                    + "  expand round: " + std::to_string(ctx. expand_round_));

        // begin v-spawn expand patten and deliver all
        // expanded patterns to other workers
        if (ctx.current_round_ > ctx.expand_round_) {
          // has reached the expand round limit
          return;
        }

        util::Output(" round: "      + std::to_string(ctx.current_round_)
                   + " generated : " + std::to_string(ctx.expand_level_.size())
                   + " pattern.");

        if (ctx.current_round_ > 0) {
          // not the initial round, all current pattern have been 
          // evaluated before, expand
          ExpandTreeLevelType next_expand_level;

          auto expand_begin = std::chrono::system_clock::now();

          // to mark where each generated pattern expanded from
          std::vector<size_t> expand_from_pattern;

          // vertical spawning
          for (int expand_level_idx = 0; 
                   expand_level_idx < static_cast<int>(ctx.expand_level_.size()); 
                   expand_level_idx++) {
            util::Debug("expanding " + std::to_string(expand_level_idx) + "'th pattern");
            util::Debug("\tsize before expand: " + std::to_string(next_expand_level.size()));

            // ###########################################################
            // ##  for the non-link pattern with no legal rhs literals  ##
            // ##  there is no way for them to expand legal gars.       ##
            // ##  if there are not too many patterns, those "hopless"  ##
            // ##  patterns would be still expanded for pruning         ##
            // ###########################################################
            bool only_expand_link = false;
            if (next_expand_level.size() > kPatternTooMany) {
              // there has already been too many patterns,
              // would not expand pattern for pruning
              if (ctx.expand_level_.node(expand_level_idx)
                                   .rhs_literals().empty()) {
                if (!GUNDAM::IsPath<true>(ctx.expand_level_.node(expand_level_idx)
                                                           .pattern())) {
                  // is not a link, won't be considered
                  continue;
                }
                // is a link but has no legal rhs literal, only need to
                // expand it as a link
                only_expand_link = true;
              }
            }
            ExpandPattern(ctx.expand_level_.node(expand_level_idx), 
                          ctx.graph_basic_statistics_, 
                          ctx. current_round_ + ctx.root_pattern_max_edge_id_,
                          ctx.restriction_,
                              only_expand_link,
                              next_expand_level);

            assert(expand_from_pattern.size() <= next_expand_level.size());
            expand_from_pattern.resize(next_expand_level.size(), expand_level_idx);
          }

          util::Output("expand end!, total pattern number: "
                    + std::to_string(next_expand_level.size()));

          auto expand_end = std::chrono::system_clock::now();

          double expand_time = CalTime(expand_begin, expand_end);

          util::Output("expand time is " + std::to_string(expand_time) + "s");

          if (next_expand_level.size() == 0) {
            /// does not find new pattern, discover end
            return;
          }
          util::Info("begin to unique pattern ");
          auto unique_begin = std::chrono::system_clock::now();
          // UniqueExpandTreeLevel(next_expand_level);
          MergeExpandedPatterns(next_expand_level,
                                ctx.expand_level_,
                                    expand_from_pattern,
                                ctx.graph_basic_statistics_,
                                ctx.ml_edge_label_set_, 
                                ctx.ml_edge_type_set_,
                                ctx.ml_to_normal_edge_label_,
                                ctx.restriction_);

          auto unique_end = std::chrono::system_clock::now();
          double unique_time = CalTime(unique_begin, unique_end);
          util::Output("unique time is " 
                    + std::to_string(unique_time) + "s");
          util::Output("uniqued pattern number: " 
                    + std::to_string(next_expand_level.size()));

          // index the uniqued patterns
          ctx.expand_level_.Swap(next_expand_level);
        }
        
        ctx.current_round_++;
        ctx.expand_level_.ReorderNodeId();

        ctx.pattern_match_count_.clear();
        ctx.pattern_match_count_.resize(ctx.expand_level_  .size(), 
                  std::vector<uint64_t>(ctx.data_graph_set_.size(), 0));

        assert(ctx.fid_ == kExpandPatternFragID);

        // ######################################################
        // ##  deliever the generated patterns to each worker  ##
        // ##  based on the number of process they have        ##
        // ######################################################
        std::vector<int> worker_id;
        for (const auto& process_num : ctx.process_num_) {
          for (int i = 0; i < process_num.second; i++)
            worker_id.emplace_back(process_num.first);
        }
        std::random_shuffle ( worker_id.begin(), 
                              worker_id.end() );

        std::vector<std::string> msg_to_workers;
        msg_to_workers.resize(ctx.frag_num_);
        // process_pattern_id_set_ is utilized as an optimization for 
        // kExpandPatternFragID to deliever message to itself
        ctx.process_pattern_id_set_.clear();
        // generate the message to deliver to workers
        for (int expand_level_idx = 0; 
                 expand_level_idx < static_cast<int>(ctx.expand_level_.size()); 
                 expand_level_idx++) {
          auto& next_expand_tree_node = ctx.expand_level_.node(expand_level_idx);
          int to_worker = worker_id.at(std::rand() % worker_id.size());
          if (next_expand_tree_node.rhs_literals().empty()) {
            // ########################################################################
            // ##  this pattern does not have rhs literal, it would only be used     ##
            // ##  to evaluate whether it is legal (support >= 1) for pruning.       ##
            // ##  if there are not many patterns, only the link would be evaluated  ##
            // ########################################################################
            if (ctx.expand_level_.size() >= kPatternTooMany) {
              if (!GUNDAM::IsPath<true>(next_expand_tree_node.pattern())) {
                // there are too many patterns, this pattern would not be evaluated
                continue;
              }
              // still need to be evaluated since it is a link even if there are 
              // already too many patterns, otherwise might miss some rules
            }
            if (to_worker == kExpandPatternFragID) {
              // as an optimization for the message deliever to itself
              //
              // "true" stands for only need to evaluate whether is legal
              ctx.process_pattern_id_set_.emplace_back(expand_level_idx, true); 
              continue;
            }
            // generate message to deliever to other workers

            // "true" stands for only need to evaluate whether is legal
            std::tuple<GraphPatternType, size_t, bool> 
              temp_message(next_expand_tree_node.pattern(), expand_level_idx, true); 

            // serialize
            msg_to_workers[to_worker] << temp_message;

            #ifndef NDEBUG
            // to test whether the serilize process is legal
            std::string temp_message_str;
            temp_message_str << temp_message;
            std::tuple<GraphPatternType, size_t, bool> temp_temp_message;
            temp_message_str >> temp_temp_message;
            std::string temp_temp_message_str;
            temp_temp_message_str << temp_temp_message;

            assert(GUNDAM::SamePattern(std::get<0>(     temp_message),
                                       std::get<0>(temp_temp_message)));
            assert(std::get<1>(     temp_message)
                == std::get<1>(temp_temp_message));
            assert(std::get<2>(     temp_message)
                == std::get<2>(temp_temp_message));
            assert(temp_message_str == temp_temp_message_str);
            #endif // NDEBUG
            continue;
          }
          // has legal rhs literal, needs to evaluate whether satisfy support bound
          if (to_worker == kExpandPatternFragID) {
            // as an optimization for the message deliever to itself
            //
            // "false" stands for it needs to evalaute whether satisfy support bound
            ctx.process_pattern_id_set_.emplace_back(expand_level_idx, false);
            continue;
          }
          // "false" stands for it needs to evalaute whether satisfy support bound
          std::tuple<GraphPatternType, size_t, bool> 
            temp_message(next_expand_tree_node.pattern(), expand_level_idx, false);
          msg_to_workers[to_worker] << temp_message;
          
          #ifndef NDEBUG
          // to test whether the serilize process is legal
          std::string temp_message_str;
          temp_message_str << temp_message;
          std::tuple<GraphPatternType, size_t, bool> temp_temp_message;
          temp_message_str >> temp_temp_message;
          std::string temp_temp_message_str;
          temp_temp_message_str << temp_temp_message;

          assert(GUNDAM::SamePattern(std::get<0>(     temp_message),
                                      std::get<0>(temp_temp_message)));
          assert(std::get<1>(     temp_message)
              == std::get<1>(temp_temp_message));
          assert(std::get<2>(     temp_message)
              == std::get<2>(temp_temp_message));
          assert(temp_message_str == temp_temp_message_str);
          #endif // NDEBUG
        }
        // send message to ordinary workers
        auto& channel_0 = messages.Channels()[0];
        for (int dst_fid = 0; dst_fid < ctx.frag_num_; dst_fid++) {
          std::string msg(kDeliverPatternPrefix);
          // would not send message to itself
          assert(dst_fid != kExpandPatternFragID
              || msg_to_workers[dst_fid] == "");
          channel_0.SendToFragment(dst_fid, std::move(msg) 
                                          + std::move(msg_to_workers[dst_fid]));
        }
        // has delieverd all message to works 
        return;
      }

      if (receive_update_pattern) {
        // ##################################################################
        // ##  receive the match number of the patterns from the ordinary  ##
        // ##  workers then generate message to ordinary worker            ##
        // ##################################################################
        ExpandTreeLevelType next_expand_level;
        // to mark whether each pattern satisfy the support bound
        std::vector<std::vector<bool>> satisfy_supp_bound;
        // the match number for each pattern on each data graph from
        // ordinary workers has already been collected in ctx.pattern_match_count_
        assert(ctx.expand_level_.size() 
            == ctx.pattern_match_count_.size());
        // #####################################################################
        // ##  first find whether it is legal (match count >= 1) on at least  ##
        // ##  one data graph                                                 ##
        // ##    if it is illegal on all data graphs, then skip this pattern  ##
        // ##    if it is legal on at least one data graph, add this pattern  ##
        // ##      to the next level of expand tree, then count the match     ##
        // ##      count for this pattern on each data graph                  ##
        // #####################################################################
        for (size_t expand_node_idx = 0; 
                    expand_node_idx < ctx.expand_level_.size(); 
                    expand_node_idx++) {
          bool is_illegal = true;
          // iterate over all data grpahs to find whether this pattern
          // is legal on at least one data graph
          for (const auto& counter : ctx.pattern_match_count_[expand_node_idx]) {
            if (counter == 0) {
              // does not legal on this data graph, try the next data graph
              continue;
            }
            // found this pattern legal on at least one data graph
            is_illegal = false;
            break;
          }
          if (is_illegal) {
            // this pattern is not legal on all data graphs
            // does not need to be considered
            continue;
          }
          // this pattern is legal on at least one data graph, add
          // to the next level of the expand tree
          next_expand_level.AddExpandTreeNode(ctx.expand_level_.node(expand_node_idx));
          // to find whether satisfy the support bound on each data graph
          auto& satisfy_supp_bound_ref
              = satisfy_supp_bound.emplace_back(std::vector<bool>());
          assert(next_expand_level.size() == satisfy_supp_bound.size());
          satisfy_supp_bound_ref.resize(ctx.pattern_match_count_[expand_node_idx].size(), false);

          for (size_t data_graph_idx = 0; 
                      data_graph_idx < ctx.pattern_match_count_[expand_node_idx].size(); 
                      data_graph_idx++) {
              
            const auto& counter = ctx.pattern_match_count_[expand_node_idx]
                                                          [ data_graph_idx];

            if (ctx.support_bound_ == 1) {
              // cannot tell the difference between whether this pattern is used to
              // only verify is legal or not or whether satisfy support bound through
              // the support counter
              if (!ctx.expand_level_.node(expand_node_idx).rhs_literals().empty()){
                // has rhs literal to verify
                satisfy_supp_bound_ref[data_graph_idx] = true;
                continue;
              }
              // does not have rhs literal to verify
              satisfy_supp_bound_ref[data_graph_idx] = false;
              continue;
            }
            if (counter == ctx.support_bound_) {
              // satisfy support bound
              satisfy_supp_bound_ref[data_graph_idx] = true;
              assert(!ctx.expand_level_.node(expand_node_idx).rhs_literals().empty());
              continue;
            }
            // does not satisfy support bound
            satisfy_supp_bound_ref[data_graph_idx] = false;
          }
          assert(satisfy_supp_bound.back().size() 
              == satisfy_supp_bound_ref.size());
          assert(satisfy_supp_bound.back().size() 
              == ctx.pattern_match_count_[expand_node_idx].size());
        }
        // whether each pattern satisfy the support bound on each
        // data graph has been generated
        assert(satisfy_supp_bound.size() 
             == next_expand_level.size());
        // index the expand_nodes
        next_expand_level.ReorderNodeId();
        ctx.expand_level_.Swap(next_expand_level);
        assert(satisfy_supp_bound.size() 
             == ctx.expand_level_.size());
        #ifndef NDEBUG
        for (int i = 0 ; i < ctx.expand_level_.size(); i++) {
          assert(ctx.expand_level_.node(i).id() == i);
        }
        #endif // NDEBUG
        // ###############################################################
        // ##  deliever the nodes in the next level of expanded tree    ##
        // ##  to each worker based on the number of process they have  ##
        // ###############################################################
        std::vector<std::string> msg_to_workers;
        msg_to_workers.resize(ctx.frag_num_);
        // an optimization, the kExpandPatternFragID no longer needs to 
        // transmit the expand node to process to itself, only needs to
        // preserve the id of those nodes
        ctx.process_expand_node_id_set_.clear();
        std::vector<std::tuple<uint32_t, // allocated rhs number
                               uint32_t, // processor(thread) number
                               uint32_t> // worker id
                   > worker_rhs_number;
        for (const auto& [worker_id, process_num] : ctx.process_num_) {
          worker_rhs_number.emplace_back(0, process_num, worker_id);
        }
        assert(worker_rhs_number.size() 
             == ctx.process_num_.size());
        for (int expand_node_idx = 0 ; 
                 expand_node_idx < ctx.expand_level_.size(); 
                 expand_node_idx++) {
          bool satisfy_supp_bound_on_at_least_one_data_graph = false;
          auto& expand_node = ctx.expand_level_.node(expand_node_idx);
          for (const bool& satisfy_supp_bound_on_data_graph
                         : satisfy_supp_bound[expand_node_idx]) {
            if (!satisfy_supp_bound_on_data_graph) {
              continue;
            }
            // satisfy support bound on at least one data graph
            satisfy_supp_bound_on_at_least_one_data_graph = true;
            break;
          }
          if (!satisfy_supp_bound_on_at_least_one_data_graph) {
            // would only be used for pruning in the expand in the next level,
            // won't be delievered to other workers for horizontal spawning
            //
            // all rhs literals on this pattern are illegal
            expand_node.rhs_literals().clear();
            continue;
          }

          assert(!expand_node.rhs_literals().empty());
          std::get<0>(worker_rhs_number.front()) += expand_node
                                                   .rhs_literals()
                                                   .size();
          
          int to_worker = std::get<2>(worker_rhs_number.front());

          std::sort(worker_rhs_number.begin(), 
                    worker_rhs_number.end(), 
                    [](const std::tuple<uint32_t, uint32_t, uint32_t>& a, 
                       const std::tuple<uint32_t, uint32_t, uint32_t>& b) -> bool { 
                    return ((double)std::get<0>(a))/((double)std::get<1>(a))
                         < ((double)std::get<0>(b))/((double)std::get<1>(b)); 
                });

          #ifndef NDEBUG 
          for (const auto& [rhs_number, 
                        process_number, 
                         worker_id] : worker_rhs_number){
            util::Debug("#" + std::to_string(rhs_number)
                      + "," + std::to_string(process_number)
                      + "," + std::to_string(worker_id)
                      + "#");
          }
          #endif // NDEBUG

          if (to_worker == kExpandPatternFragID) {
            // to itself, does not need message transformation
            ctx.process_expand_node_id_set_.emplace_back(expand_node_idx);
            continue;
          }
          // to other workers, generated message
          util::Debug("processing expand node: " + std::to_string(expand_node.id())
                                + " to worker: " + std::to_string(to_worker));
          msg_to_workers[to_worker] << expand_node;
        }
        auto& channel_0 = messages.Channels()[0];
        for (int dst_fid = 0; dst_fid < ctx.frag_num_; dst_fid++) {
          std::string msg(kDeliverExpandNodesPrefix);
          channel_0.SendToFragment(dst_fid, std::move(msg) 
                                          + std::move(msg_to_workers[dst_fid]));
        }
        return;
      }
    }
    else {
      assert(receive_deliver_expand_node
          || receive_deliver_pattern);
      // oridinary node receive deliver message
      if (receive_deliver_expand_node)  {
        ctx.process_expand_node_id_set_.clear();
        ctx.process_expand_node_id_set_.reserve(ctx.expand_level_.size());
        for (int i = 0; i < ctx.expand_level_.size(); i++){
          ctx.process_expand_node_id_set_.emplace_back(i);
        }
      }
    }

    // ordinary workers should not receive update_expand_node or update_pattern
    assert(receive_deliver_expand_node
        || receive_deliver_pattern);

    if (receive_deliver_pattern) {
      // ###################################################################
      // ##  receive the pattern delivered from the kExpandPatternFragID  ##
      // ##  calculate the match count for each pattern based on the      ##
      // ##  requirement, send the calculated match count back            ##
      // ###################################################################

      std::vector<std::tuple<size_t, // pattern_id
                           uint32_t, // pattern_supp
                             size_t> // graph_id
                 > patterns_support;

      size_t total_patterns_to_evalute = (ctx.fid_ == kExpandPatternFragID)?
                                          // for node kExpandPatternFragID, 
                                          // the id of patterns it needs to evaluate
                                          // are stored in kExpandPatternFragID
                                          ctx.process_pattern_id_set_.size()
                                          // for ordinary nodes, the patterns
                                          // for it to be evaluted are stored in
                                          // delivered_patterns
                                                 : delivered_patterns.size();

      util::Info("total patterns to evalute: "
                + std::to_string(total_patterns_to_evalute));

      omp_lock_t pattern_support_lock;
      omp_init_lock(&pattern_support_lock);
                
      omp_set_nested(false);
      #pragma omp parallel for schedule(dynamic) 
      for (int i = 0; i < total_patterns_to_evalute; i++) {
        util::Debug("* Evaluate match count: number of threads in the team - "
                  + std::to_string(omp_get_num_threads()) + "\n");

        auto kDeliveredPatternId = (ctx.fid_ == kExpandPatternFragID)?
                                    ctx.process_pattern_id_set_[i].first
                                  : std::get<1>(delivered_patterns[i]);

        util::Debug("kDeliveredPatternId: "
                   + std::to_string(kDeliveredPatternId));

        // the pattern for this worker to match
        auto& delivered_pattern = (ctx.fid_ == kExpandPatternFragID)?
                                   ctx.expand_level_.node(kDeliveredPatternId).pattern()
                                 : std::get<0>(delivered_patterns[i]);

        // whether to only verify whether it is legal
        const bool  kDeliveredPatternOnlyVerifyIsLegal 
                                 = (ctx.fid_ == kExpandPatternFragID)?
                                    ctx.process_pattern_id_set_[i].second
                                  : std::get<2>(delivered_patterns[i]);

        util::Debug("kDeliveredPatternOnlyVerifyIsLegal: "
                   + std::to_string(kDeliveredPatternOnlyVerifyIsLegal));

        #ifndef NDEBUG
        if (ctx.fid_ == kExpandPatternFragID) {
          // if only verify whether it is legal, then it should not have legal rhs literal
          assert(( kDeliveredPatternOnlyVerifyIsLegal &&  ctx.expand_level_.node(kDeliveredPatternId).rhs_literals().empty())
              || (!kDeliveredPatternOnlyVerifyIsLegal && !ctx.expand_level_.node(kDeliveredPatternId).rhs_literals().empty()));
        }
        #endif // NDEBUG

        // if only needs to verify whether it is legal, then the kMatchCountLimit
        // would be set to 1
        const size_t kMatchCountLimit = kDeliveredPatternOnlyVerifyIsLegal ?
                                        1 : ctx.support_bound_;

        // verify this pattern on all data graphs
        for (size_t data_graph_idx = 0; 
                    data_graph_idx < ctx.data_graph_set_.size(); 
                    data_graph_idx++) {

          auto& data_graph = ctx.data_graph_set_[data_graph_idx]
                                .data_graph();

          int match_count = 0;
          
          if ( ctx.restriction_.gar() 
            || ctx.restriction_.horn_rule() ) {
            // ##############################################
            // ##  both gar and horn_rule use isomorphism  ##
            // ##############################################
            match_count = GUNDAM::MatchUsingMatch(
                        delivered_pattern,
                        data_graph, 
                        kMatchCountLimit, 
                        ctx.time_limit_);
          }
          else {
            // ##########################
            // ##  gcr use Simulation  ##
            // ##########################
            std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
            std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> > candidate_set;
            assert(ctx.restriction_.gcr());
            const bool kInitCandidateSetSucc 
              = GUNDAM::_dp_iso_using_match
                      ::InitCandidateSet<GUNDAM::MatchSemantics::kHomomorphism>(
                        delivered_pattern,
                        data_graph, 
                        candidate_set);
            if (kInitCandidateSetSucc) {
              if (ctx.restriction_.central_to_central()) {
                std::vector<GraphPatternType> 
                       connected_components(GUNDAM::ConnectedComponent(delivered_pattern));
                assert(connected_components.size() == 2);
                for (const auto& connected_component
                               : connected_components) {
                  const auto [end_vertex_handle_set,
                          central_vertex_handle] = GUNDAM::StarEndPoints<true>(connected_component);
                  assert(end_vertex_handle_set.size() >= 2);
                  assert(end_vertex_handle_set.size() == 2 || central_vertex_handle);
                  std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type> central_vertex_handle_set;
                  if (end_vertex_handle_set.size() == 2) {
                    // is link, all vertexes can be legal central vertex, add them all
                    assert(GUNDAM::IsPath<true>(connected_component));
                    for (auto vertex_it = connected_component.VertexBegin();
                             !vertex_it.IsDone();
                              vertex_it++) {
                      assert(delivered_pattern.FindVertex(vertex_it->id()));
                      central_vertex_handle_set.emplace_back(delivered_pattern.FindVertex(vertex_it->id()));
                    }
                    assert(central_vertex_handle_set.size() == connected_component.CountVertex());
                  }
                  else {
                    assert(delivered_pattern.FindVertex(central_vertex_handle->id()));
                    central_vertex_handle_set.emplace_back(delivered_pattern.FindVertex(central_vertex_handle->id()));
                    assert(central_vertex_handle_set.size() == 1);
                  }
                  for (const auto& vertex_handle
                         : central_vertex_handle_set) {
                    auto central_candidate_set_it = candidate_set.find(vertex_handle);
                    assert(central_candidate_set_it
                                != candidate_set.end());
                    
                    std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> selected_central_candidate_set;
                    selected_central_candidate_set.reserve(
                             central_candidate_set_it->second.size());
                    assert(!ctx.data_graph_set_[data_graph_idx]
                                          .centrel_vertex_set(vertex_handle->label()).empty());
                    const auto& data_graph_centrel_vertex_set
                          = ctx.data_graph_set_[data_graph_idx]
                                          .centrel_vertex_set(vertex_handle->label());
                    std::set_intersection(central_candidate_set_it->second.begin(),
                                          central_candidate_set_it->second. end (),
                                             data_graph_centrel_vertex_set.begin(),
                                             data_graph_centrel_vertex_set. end (),
                                          std::back_inserter(
                                               selected_central_candidate_set));
                    selected_central_candidate_set.swap(central_candidate_set_it->second);
                  }
                }
              }
              const bool kRefineCandidateSetSucc 
                = GUNDAM::_dp_iso_using_match
                        ::RefineCandidateSet( delivered_pattern,
                                              data_graph, 
                                              candidate_set );
              if (kRefineCandidateSetSucc) {
                std::vector<GraphPatternType> connected_components
                 = GUNDAM::ConnectedComponent(delivered_pattern);
                assert(connected_components.size() == 2);
                // ##########################################################
                // ##  only need to multiply the minimal size              ##
                // ##  of two candidate sets, could know whether there     ##
                // ##  could have literals that satisfy the support bound  ##
                // ##########################################################
                std::vector<size_t> minimal_candidate_size(2, 0);
                assert(minimal_candidate_size.size() == 2);
                for (size_t cc_idx = 0;
                            cc_idx < connected_components.size();
                            cc_idx++) {
                  const auto& cc  =  connected_components[cc_idx];
                  for (auto vertex_it = cc.VertexBegin();
                           !vertex_it.IsDone();
                            vertex_it++) {
                    assert(delivered_pattern.FindVertex(vertex_it->id()));
                    assert(candidate_set.find(delivered_pattern.FindVertex(vertex_it->id()))
                        != candidate_set.end());
                    if (minimal_candidate_size[cc_idx] == 0) {
                      minimal_candidate_size[cc_idx]
                            = candidate_set [delivered_pattern.FindVertex(vertex_it->id())].size();
                      continue;
                    }
                    minimal_candidate_size[cc_idx] = minimal_candidate_size[cc_idx] < candidate_set[delivered_pattern.FindVertex(vertex_it->id())].size()?
                                                     minimal_candidate_size[cc_idx] : candidate_set[delivered_pattern.FindVertex(vertex_it->id())].size();
                  }
                }
                assert(minimal_candidate_size[0] > 0);
                assert(minimal_candidate_size[1] > 0);
                match_count = minimal_candidate_size[0] * minimal_candidate_size[1];
                match_count = match_count < kMatchCountLimit?
                              match_count : kMatchCountLimit;
              }
              #ifndef NDEBUG
              else {
                assert(match_count == 0);
              }
              #endif // NDEBUG
            }
            #ifndef NDEBUG
            else {
              assert(match_count == 0);
            }
            #endif // NDEBUG
          }
                    
          util::Debug("< kMatchCountLimit: " + std::to_string(kMatchCountLimit)
                    + ", match_count: "      + std::to_string(match_count) + " >");

          assert(match_count >= 0 && match_count <= kMatchCountLimit);

          #ifndef NDEBUG
          if (match_count == ctx.support_bound_) {
            assert(!ctx.expand_level_.node(kDeliveredPatternId).rhs_literals().empty());
          }
          #endif // NDEBUG

          omp_set_lock(&pattern_support_lock);
          patterns_support.emplace_back(kDeliveredPatternId, match_count, data_graph_idx);
          omp_unset_lock(&pattern_support_lock);
        }
      }

      if (ctx.fid_ == kExpandPatternFragID) {
        // an optimization, simply store the match count in 
        // ctx.pattern_match_count_ does not need to transform
        // the message to it self
        for (const auto& [pattern_id, 
                          pattern_supp, 
                          graph_id] : patterns_support) {
          assert(pattern_id >= 0 
              && pattern_id <  ctx.pattern_match_count_.size());
          assert(graph_id >= 0 
              && graph_id <  ctx.pattern_match_count_[pattern_id].size());
          ctx.pattern_match_count_[pattern_id][graph_id] = pattern_supp;
          #ifndef NDEBUG
          if (ctx.support_bound_ != 1
           && ctx.support_bound_ == pattern_supp) {
            assert(!ctx.expand_level_.node(pattern_id).rhs_literals().empty());
          }
          #endif // NDEBUG
        }
        auto& channel_0 = messages.Channels()[0];
        std::string msg(kUpdatePatternPrefix);
        channel_0.SendToFragment(kExpandPatternFragID, msg);
        return;
      }

      // is not kExpandPatternFragID, needs to genrate message to kExpandPatternFragID
      assert(ctx.fid_ != kExpandPatternFragID);

      auto& channel_0 = messages.Channels()[0];
      // send match count for each pattern on each data graph, back to kExpandPatternFragID
      std::string msg(kUpdatePatternPrefix);
      // the kExpandPatternFragID does not need to transform data through the
      // MPI method through itself, 
      for (const auto& [pattern_id, 
                        pattern_supp, 
                          graph_id] : patterns_support) {
        util::Debug(" pattern_id: " + std::to_string(pattern_id)
                + " pattern_supp: " + std::to_string(pattern_supp)
                    + " graph_id: " + std::to_string(graph_id));
        if (pattern_supp == 0) {
          // there is no need to pass such info
          continue;
        }
        msg = std::move(msg) + " " + std::to_string(pattern_id)
                             + " " + std::to_string(pattern_supp)
                             + " " + std::to_string(graph_id);
      }
      channel_0.SendToFragment(kExpandPatternFragID, msg);
      return;
    }

    // ordinary workers should not receive update_expand_node or update_pattern
    assert(receive_deliver_expand_node);
    
    // #######################################################################
    // ##  process the expand node delivered from the kExpandPatternFragID  ##
    // ##  expand lhs literals on each pattern                              ##
    // #######################################################################

    if (ctx.process_expand_node_id_set_.size() == 0) {
      // ##########################################
      // ## does not have expand node to process ##
      // ##########################################
      auto& channel_0 = messages.Channels()[0];
      std::string msg(kUpdateExpandNodesPrefix);
      channel_0.SendToFragment(kExpandPatternFragID, msg);
      return;
    }
    // #################################
    // ## have expand node to process ##
    // #################################    
    assert(ctx.process_expand_node_id_set_.size() > 0);
    // each level has same edge number
    const SubstructureLevelType kCurrentPatternEdgeSize
           = ctx.expand_level_.node(ctx.process_expand_node_id_set_.front())
                              .const_pattern()
                              .CountEdge();

    // current level of the generation tree,
    GenerateTreeLevelType current_generate_level;

    omp_lock_t current_generate_level_lock;
    omp_init_lock(&current_generate_level_lock);

    util::Info("number of nodes for horizontal spawning: "
              + std::to_string( ctx.process_expand_node_id_set_.size() ));

    const int kTenPercent = ctx.process_expand_node_id_set_.size() / 10.0;
    // horizontal spawning
    omp_set_nested(true);
    // omp_set_max_active_levels(2);
    #pragma omp parallel for schedule(dynamic) 
    for (int process_expand_node_idx = 0; 
             process_expand_node_idx < static_cast<int>(ctx.process_expand_node_id_set_.size()); 
             process_expand_node_idx++) {
      // add the generated pattern to the current_generate_level
      if ((process_expand_node_idx != 0) && (kTenPercent != 0) 
       && (process_expand_node_idx % kTenPercent == 0)) {
        util::Info("horizontal spawning completed: "
                 + std::to_string( process_expand_node_idx / kTenPercent )
                 + "0% ");
      }
      util::Debug("# Evaluate literal support: number of threads in the team - "
                 + std::to_string( omp_get_num_threads() )
                 + "\n");
      // get the reference of the current processing expand node
      auto& current_expand_node = ctx.expand_level_.node(
                                  ctx.process_expand_node_id_set_[
                                      process_expand_node_idx]);
      const auto kCurrentExpandNodeId
                = current_expand_node.id();

      assert(kCurrentExpandNodeId >= process_expand_node_idx);

      // get the constant reference of the pattern of the current
      // processing expand node for further processing
      const auto& current_expand_pattern
                = current_expand_node.const_pattern();

      assert(!current_expand_node.rhs_literals().empty());

      // node on the generate tree corresponds to the current
      // processing expand node
      GenerateTreeNodeType current_generate_node(kCurrentExpandNodeId,
                                                  current_expand_pattern);

      using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
      using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;
      
      using CandidateSetContainer 
        = std::map<GraphPatternVertexHandleType,
          std::vector<DataGraphVertexHandleType>>;

      // candidate set on each data graph
      std::vector<CandidateSetContainer> candidate_set_for_data_graph_set;
      candidate_set_for_data_graph_set.reserve(ctx.data_graph_set_.size());

      // candidate set removed nec on each data graph
      std::vector<CandidateSetContainer> candidate_set_removed_nec_for_data_graph_set;
      candidate_set_removed_nec_for_data_graph_set.reserve(ctx.data_graph_set_.size());

      for (size_t data_graph_idx = 0; 
                  data_graph_idx < ctx.data_graph_set_.size(); 
                  data_graph_idx++) {

        CandidateSetContainer candidate_set;

        bool init_candidate_set_succ = false;
        if (ctx.restriction_.gcr()) {
          // gcr, use homomorphism
          init_candidate_set_succ
          = GUNDAM::_dp_iso_using_match
                  ::InitCandidateSet<GUNDAM::MatchSemantics::kHomomorphism>(
            current_generate_node.pattern(),
            ctx.data_graph_set_[data_graph_idx].data_graph(), 
            candidate_set);
        }
        else {
          // gar or horn rule, use isomorphism
          init_candidate_set_succ
          = GUNDAM::_dp_iso_using_match
                  ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
            current_generate_node.pattern(),
            ctx.data_graph_set_[data_graph_idx].data_graph(), 
            candidate_set);
        }

        if (!init_candidate_set_succ) {
          // assert(false); can be empty on one of these data graphs
          //
          // emplace back an empty candidate set
          candidate_set_for_data_graph_set.emplace_back();
          candidate_set_removed_nec_for_data_graph_set.emplace_back();
          continue;
        }
        
        if (ctx.restriction_.central_to_central()) {
          std::vector<GraphPatternType> 
                connected_components(GUNDAM::ConnectedComponent(current_generate_node.pattern()));
          assert(connected_components.size() == 2);
          for (const auto& connected_component
                         : connected_components) {
            const auto [end_vertex_handle_set,
                    central_vertex_handle] = GUNDAM::StarEndPoints<true>(connected_component);
            assert(end_vertex_handle_set.size() >= 2);
            assert(end_vertex_handle_set.size() == 2 || central_vertex_handle);
            std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type> central_vertex_handle_set;
            if (end_vertex_handle_set.size() == 2) {
              // is link, all vertexes can be legal central vertex, add them all
              assert(GUNDAM::IsPath<true>(connected_component));
              for (auto vertex_it = connected_component.VertexBegin();
                        !vertex_it.IsDone();
                        vertex_it++) {
                assert(current_generate_node.pattern().FindVertex(vertex_it->id()));
                central_vertex_handle_set.emplace_back(current_generate_node.pattern().FindVertex(vertex_it->id()));
              }
              assert(central_vertex_handle_set.size() == connected_component.CountVertex());
            }
            else {
              assert(current_generate_node.pattern().FindVertex(central_vertex_handle->id()));
              central_vertex_handle_set.emplace_back(current_generate_node.pattern().FindVertex(central_vertex_handle->id()));
              assert(central_vertex_handle_set.size() == 1);
            }
            for (const auto& vertex_handle
                   : central_vertex_handle_set) {
              auto central_candidate_set_it = candidate_set.find(vertex_handle);
              assert(central_candidate_set_it
                          != candidate_set.end());
              
              std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> selected_central_candidate_set;
              selected_central_candidate_set.reserve(
                        central_candidate_set_it->second.size());
              assert(!ctx.data_graph_set_[data_graph_idx]
                                    .centrel_vertex_set(vertex_handle->label()).empty());
              const auto& data_graph_centrel_vertex_set
                    = ctx.data_graph_set_[data_graph_idx]
                                    .centrel_vertex_set(vertex_handle->label());
              std::set_intersection(central_candidate_set_it->second.begin(),
                                    central_candidate_set_it->second. end (),
                                       data_graph_centrel_vertex_set.begin(),
                                       data_graph_centrel_vertex_set. end (),
                                    std::back_inserter(
                                         selected_central_candidate_set));
              selected_central_candidate_set.swap(central_candidate_set_it->second);
            }
          }
        }

        const bool kRefineCandidateSetSucc 
          = GUNDAM::_dp_iso_using_match
                  ::RefineCandidateSet(
            current_generate_node.pattern(),
            ctx.data_graph_set_[data_graph_idx].data_graph(), 
            candidate_set);

        if (!kRefineCandidateSetSucc) {
          // assert(false); can be empty on one of these data graphs
          //
          // emplace back an empty candidate set
          candidate_set_for_data_graph_set.emplace_back();
          candidate_set_removed_nec_for_data_graph_set.emplace_back();
          continue;
        }
        // add the refined candidate set to candidate_set_for_data_graph_set
        candidate_set_for_data_graph_set.emplace_back(candidate_set); // does not use move, since candidate_set would be used later

        GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                            DataGraphType>(
                  candidate_set,
                  ctx.data_graph_set_[data_graph_idx]
                     .data_graph_nec());
        // if it is not empty before remove Nec, it should not be empty after remove it
        assert(!candidate_set.empty());
        
        candidate_set_removed_nec_for_data_graph_set.emplace_back(std::move(candidate_set));
      }
      assert(            candidate_set_for_data_graph_set.size() == ctx.data_graph_set_.size());
      assert(candidate_set_removed_nec_for_data_graph_set.size() == ctx.data_graph_set_.size());
      
      // all existed rhs literal have been added into rhs_literal_info_set
      // clear the current_expand_node.rhs_literals()
      std::vector<DataGraphLiteralType> rhs_literal_info_set;
      rhs_literal_info_set.swap(current_expand_node.rhs_literals());
      util::Debug("rhs_literal_info_set.size(): "
                 + std::to_string(rhs_literal_info_set.size()));
      std::sort(rhs_literal_info_set.begin(),
                rhs_literal_info_set.end());
      assert(std::is_sorted(rhs_literal_info_set.begin(),
                            rhs_literal_info_set.end()));

      #ifndef NDEBUG
      for (const auto& literal : rhs_literal_info_set) {
        assert(literal.literal_type() != gar::LiteralType::kMlLiteral);
      }
      #endif // NDEBUG

      std::vector<DataGraphLiteralType> lhs_literal_info_not_in_rhs_set;

      assert(ctx.j_ > 1 || current_expand_node.lhs_literals().empty());
      if (ctx.j_ > 1) {


        lhs_literal_info_not_in_rhs_set.swap(current_expand_node.lhs_literals());

        util::Debug("lhs literal size before expand: "
                   + std::to_string(lhs_literal_info_not_in_rhs_set.size()));

        ExpandLiteral(current_expand_node, 
                      ctx.graph_basic_statistics_,
                      ctx.ml_edge_label_set_,
                      ctx.ml_edge_type_set_,
                      ctx.ml_to_normal_edge_label_,
                      ctx.restriction_,
                      lhs_literal_info_not_in_rhs_set);
        
        util::Debug("lhs literal size after expand: "
                   + std::to_string(lhs_literal_info_not_in_rhs_set.size()));

        std::sort(lhs_literal_info_not_in_rhs_set.begin(),
                  lhs_literal_info_not_in_rhs_set.end());
        assert(std::is_sorted(lhs_literal_info_not_in_rhs_set.begin(),
                              lhs_literal_info_not_in_rhs_set.end()));
        
        lhs_literal_info_not_in_rhs_set.erase(
            std::unique( lhs_literal_info_not_in_rhs_set.begin(), 
                         lhs_literal_info_not_in_rhs_set.end() ), 
                         lhs_literal_info_not_in_rhs_set.end() );
        
        util::Debug("lhs literal size after unique: "
                   + std::to_string(lhs_literal_info_not_in_rhs_set.size()));

        for (auto lhs_literal_it  = lhs_literal_info_not_in_rhs_set.begin();
                  lhs_literal_it != lhs_literal_info_not_in_rhs_set.end();) {
          if (lhs_literal_it->literal_type() == gar::LiteralType::kEdgeLiteral
          || !ctx.restriction_.LiteralTypeConsideredInLhs(lhs_literal_it->literal_type())) {
            // does not add edge literal in lhs
            lhs_literal_it = lhs_literal_info_not_in_rhs_set.erase(lhs_literal_it);
            continue;
          }
          if (std::binary_search(rhs_literal_info_set.begin(),
                                 rhs_literal_info_set.end(),
                                *lhs_literal_it)) {
            // has_same_in_rhs
            lhs_literal_it = lhs_literal_info_not_in_rhs_set.erase(lhs_literal_it);
            continue;
          }
          // preserve
          lhs_literal_it++;
        }
        
        util::Debug("lhs literal size after remove literal in rhs: "
                   + std::to_string(lhs_literal_info_not_in_rhs_set.size()));
      }

      // they should be empty now
      assert(current_expand_node.RhsLiteralsBegin()
          == current_expand_node.RhsLiteralsEnd());
      assert(current_expand_node.LhsLiteralsBegin()
          == current_expand_node.LhsLiteralsEnd());

      // first to find whether that are literals automorphism to each other
      // use disjoint_set to merge literals that are automorphism to each other
      GUNDAM::DisjointSet<size_t> rhs_disjoint_set(rhs_literal_info_set.size());
      std::vector<std::pair<GarType, size_t>> existed_gars;
      for (size_t rhs_literal_idx = 0; 
                  rhs_literal_idx < rhs_literal_info_set.size(); 
                  rhs_literal_idx++) {
        const auto& literal_info = rhs_literal_info_set[rhs_literal_idx];
        GarType new_gar(current_generate_node.const_pattern());
        new_gar.AddY(literal_info);
        bool has_automorphism = false;
        for (auto& [existed_gar,
                    existed_gar_rhs_literal_idx]
                  : existed_gars) {
          if (!gar::SameGar(existed_gar, new_gar)) {
            continue;
          }
          has_automorphism = true;
          rhs_disjoint_set.Merge(rhs_literal_idx,
                     existed_gar_rhs_literal_idx);
          break;
        }
        if (has_automorphism) {
          continue;
        }
        // does not have automorphism literals
        existed_gars.emplace_back(new_gar, rhs_literal_idx);
      }
      existed_gars.clear();

      // rhs_disjoint_set holds all rhs literals that are automorphism
      // to each other, select one representative literal from each set
      std::vector<size_t> rhs_literal_id_set_to_evaluate;
      for (size_t rhs_literal_info_idx = 0; 
                  rhs_literal_info_idx < rhs_literal_info_set.size(); 
                  rhs_literal_info_idx++) {
        size_t group_id = rhs_disjoint_set.Find(rhs_literal_info_idx);
        rhs_literal_id_set_to_evaluate.emplace_back(group_id);
      }
      util::Debug("rhs literal number before remove automorphism: " 
                 + std::to_string(rhs_literal_info_set.size()));
      std::sort( rhs_literal_id_set_to_evaluate.begin(), 
                 rhs_literal_id_set_to_evaluate.end() );
      rhs_literal_id_set_to_evaluate.erase(
          std::unique( rhs_literal_id_set_to_evaluate.begin(), 
                       rhs_literal_id_set_to_evaluate.end() ), 
                       rhs_literal_id_set_to_evaluate.end() );
      util::Debug("rhs literal number after remove automorphism: " 
                 + std::to_string(rhs_literal_id_set_to_evaluate.size()));

      std::vector<bool> 
         literal_can_be_considered_in_rhs(rhs_literal_info_set.size(), false), // staisfy support bound
         literal_can_be_considered_in_lhs(rhs_literal_info_set.size(), false), // legal, supp >= 1
         // to mark whether this representive literal have been considered
                   rhs_literal_considered(rhs_literal_info_set.size(), false);

      std::vector<float> rhs_literal_conf(rhs_literal_info_set.size(), -1.0);
      std::vector<int> rhs_literal_supp(rhs_literal_info_set.size(), 0);
      std::vector<double> rhs_literal_probability(rhs_literal_info_set.size(), 0.0);
      std::vector<std::string> root_gars_to_name(rhs_literal_info_set.size());
      uint64_t gar_offset = 0;

      for (size_t rhs_literal_id_set_to_evaluate_idx = 0; 
                  rhs_literal_id_set_to_evaluate_idx < rhs_literal_id_set_to_evaluate.size(); 
                  rhs_literal_id_set_to_evaluate_idx++) {
        size_t group_id = rhs_literal_id_set_to_evaluate[rhs_literal_id_set_to_evaluate_idx];
        assert(group_id < rhs_literal_info_set.size());
        assert(group_id >= 0);
        if (rhs_literal_considered[group_id]) {
          continue;
        }
        rhs_literal_considered[group_id] = true;

        const auto& literal_info = rhs_literal_info_set[group_id];

        // #########################################
        // ##  hold the set of literals as X U Y  ##
        // #########################################
        const std::vector<DataGraphLiteralType> literal_set = {literal_info};

        assert(candidate_set_for_data_graph_set .size()
                          == ctx.data_graph_set_.size());

        assert(candidate_set_removed_nec_for_data_graph_set .size()
                                      == ctx.data_graph_set_.size());

        uint64_t supp = 0;
          float  conf = 0;
        double probability = 0.0;
        std::vector<std::pair<std::vector<VertexIDType>, std::vector<EdgeIDType>>> match_vec;
        // ######################################################################
        // ##  evaluate gar on all data graphs until it satisfy support bound  ##
        // ##  on one of these data graphs                                     ##
        // ######################################################################
        for (size_t data_graph_idx = 0; 
                    data_graph_idx < ctx.data_graph_set_.size(); 
                    data_graph_idx++) {

          assert(data_graph_idx < ctx.data_graph_set_.size());
          auto& data_graph = ctx.data_graph_set_.at(data_graph_idx).data_graph();

          assert(data_graph_idx < candidate_set_for_data_graph_set.size());
          const auto& candidate_set
                    = candidate_set_for_data_graph_set.at(data_graph_idx);
          if (candidate_set.empty()) {
            assert(candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx).empty());
            continue;
          }

          assert(data_graph_idx < candidate_set_removed_nec_for_data_graph_set.size());
          const auto& candidate_set_removed_nec
                    = candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx);
          if (candidate_set_removed_nec.empty()) {
            assert(false);
            continue;
          }

          auto [temp_supp, temp_conf] 
             = HasMatch(current_generate_node.pattern(),
                        literal_set,
                        literal_info,
                        data_graph,
                        ctx.normal_to_ml_edge_label_,
                        candidate_set,
                        candidate_set_removed_nec,
                        ctx.support_bound_, false,
                        ctx.graph_basic_statistics_,
                        ctx.restriction_,
                        ctx.time_limit_,
                        ctx.time_limit_per_supp_,
                        ctx.store_match_,
                        match_vec,
                        probability);
          //std::cout << "line 2274 prob " << probability << std::endl;
          std::string literal_str;
          literal_str << literal_info;

          std::string pattern_str;
          pattern_str << current_generate_node.pattern();

          util::Debug( "literal_str: " + literal_str );
          util::Debug( "pattern_str: " + pattern_str );
          util::Debug( "temp_supp: " + std::to_string(temp_supp) );
          util::Debug( "temp_conf: " + std::to_string(temp_conf) );

          assert(temp_supp >= 0
             && (temp_supp <= ctx.support_bound_
              || ctx.restriction_.specified_confidence_bound()));

          if (temp_supp == 0) {
            // is not legal on this data graph
            continue;
          }
          if (supp == 0) {
            supp = temp_supp;
          }

          assert(ctx.restriction_.specified_confidence_bound()
              || supp <= ctx.support_bound_);
          if (!ctx.restriction_.specified_confidence_bound()) {
            if (supp == ctx.support_bound_) {
              // statisfy support bound on this data graph
              // does not need to consider other data graphs
              break;
            }
            continue;
          }
          if (temp_supp < ctx.support_bound_) {
            continue;
          }
          supp = temp_supp;
          if (temp_conf >= ctx.restriction_.confidence_bound()) {
            conf = temp_conf;
            break;
          }
        }
        assert(conf >= ctx.restriction_.confidence_bound()
            || conf == 0);

        // this literal has not been considered before
        assert(!literal_can_be_considered_in_rhs[group_id]);
        assert(!literal_can_be_considered_in_lhs[group_id]);

        if (supp > 0) {
          // is legal on at least one data graph
          // but not satisfy support bound on all data graphs
          // not satisfy support bound but is legal
          if (literal_info.literal_type() != gar::LiteralType::kEdgeLiteral) {
            literal_can_be_considered_in_lhs[group_id] = true;
          }
        }

        if (supp >= ctx.support_bound_) {
          // has statisfied support bound on at least one data graph
          literal_can_be_considered_in_rhs[group_id] = true;
          rhs_literal_supp[group_id] = supp;
          rhs_literal_probability[group_id] = probability;
          if ( ctx.restriction_.specified_confidence_bound() ) {
            rhs_literal_conf[group_id] = conf;
          }
          if (ctx.store_match_) {
            std::string root_gar_name = ExportGARMatch(literal_info, match_vec, kCurrentPatternEdgeSize,
                        process_expand_node_idx, gar_offset, ctx.fid_, ctx.output_match_dir_);
            gar_offset++;
            root_gars_to_name[group_id] = root_gar_name;
          }
        }
      }

      assert(current_expand_node.RhsLiteralsBegin()
          == current_expand_node.RhsLiteralsEnd());

      // to hold the literals that are not the representative one 
      // in the group, their support would not be calculated
      std::vector<DataGraphLiteralType> 
        rhs_literal_info_satisfy_supp_bound_set_duplicated;

      std::vector<DataGraphLiteralType> 
        current_expand_node_rhs_literal;
      current_expand_node_rhs_literal.reserve(
                          rhs_literal_info_set.size());

      std::vector<float> current_expand_node_rhs_literal_conf;
      std::vector<int> current_expand_node_rhs_literal_supp;
      std::vector<std::string> current_expand_node_rhs_literal_name;
      std::vector<double> current_expand_node_rhs_literal_probability;
      if (ctx.restriction_.specified_confidence_bound()) {
        current_expand_node_rhs_literal_conf.reserve(
                            rhs_literal_info_set.size());
      }
      current_expand_node_rhs_literal_supp.reserve(
                          rhs_literal_info_set.size());
      current_expand_node_rhs_literal_probability.reserve(
                          rhs_literal_info_set.size());
      if (ctx.store_match_) {
        current_expand_node_rhs_literal_name.reserve(
                          rhs_literal_info_set.size());
      }

      // for the literals legal satisfy support bound
      //    add it to both rhs and lhs literals set, 
      // for the literals legal on data graph but not satisfy support bound
      //    add it to lhs literals set only
      for (size_t rhs_literal_info_idx = 0; 
                  rhs_literal_info_idx < rhs_literal_info_set.size(); 
                  rhs_literal_info_idx++) {
        size_t group_id = rhs_disjoint_set.Find(rhs_literal_info_idx);
        const auto& rhs_literal_info
                  = rhs_literal_info_set[rhs_literal_info_idx];
        if (literal_can_be_considered_in_rhs[group_id]) {
          assert(rhs_literal_info.literal_type() == gar::LiteralType::kEdgeLiteral 
              || literal_can_be_considered_in_lhs[group_id]);
          if (ctx.j_ > 1) {
            // is legal, add to Lhs
            if (rhs_literal_info.literal_type() != gar::LiteralType::kEdgeLiteral) {
              // edge literal won't be added into lhs
              if (ctx.restriction_.LiteralTypeConsideredInLhs(rhs_literal_info.literal_type())) {
                current_expand_node.AddLhsLiteral(rhs_literal_info);
              }
            }
          }
          // satisfy the support bound
          if (group_id == rhs_literal_info_idx) {
            // is the representive element of this group, add to rhs
            // would be considerd
            if (ctx.restriction_.specified_confidence_bound()) {
              assert(rhs_literal_conf[group_id] >= 0);
              current_expand_node_rhs_literal_conf.emplace_back(
                                  rhs_literal_conf[group_id]);
            }
            current_expand_node_rhs_literal_supp.emplace_back(
                                  rhs_literal_supp[group_id]);
            current_expand_node_rhs_literal_probability.emplace_back(
                                  rhs_literal_probability[group_id]);
            if (ctx.store_match_) {
              current_expand_node_rhs_literal_name.emplace_back(
                                 root_gars_to_name[group_id]);
            }
            current_expand_node_rhs_literal.emplace_back(rhs_literal_info);
            continue;
          }
          // is not the representive element of this group
          // add to rhs_literal_info_satisfy_supp_bound_set_duplicated
          // won't be used to initialize literal trees, would be added
          // back to expand node since they might nolonger automorphism
          // to each other after expandation
          if (ctx.restriction_.specified_confidence_bound()
           && ctx.restriction_.confidence_bound() <= rhs_literal_conf[group_id]) {
            // such a literal does not satisfy the confidence bound, does not
            // need to be considered in rhs_literal_info_satisfy_supp_bound_set_duplicated
            continue;
          }
          rhs_literal_info_satisfy_supp_bound_set_duplicated.emplace_back(
                              rhs_literal_info);
          continue;
        }
        if (ctx.j_ == 1) {
          // does not need to add lhs literals
          continue;
        }
        if (literal_can_be_considered_in_lhs[group_id]) {
          // does not satisfy the support bound but is legal,
          // add to Lhs literal set
          if (rhs_literal_info.literal_type() != gar::LiteralType::kEdgeLiteral) {
            // edge literal won't be added into lhs
            if (ctx.restriction_.LiteralTypeConsideredInLhs(rhs_literal_info.literal_type())) {
              current_expand_node.AddLhsLiteral(rhs_literal_info);
            }
          }
        }
      }

      if (ctx.j_ > 1) {
        // consider lhs literal only when j > 1        
        GUNDAM::DisjointSet<size_t> lhs_disjoint_set(
                                    lhs_literal_info_not_in_rhs_set.size());
        std::vector<std::pair<GarType, size_t>> existed_gars;
        for (size_t lhs_literal_idx = 0; 
                    lhs_literal_idx < lhs_literal_info_not_in_rhs_set.size(); 
                    lhs_literal_idx++) {
          const auto& literal_info
                = lhs_literal_info_not_in_rhs_set[lhs_literal_idx];
          const std::vector<DataGraphLiteralType> literal_set = {literal_info};
          GarType new_gar(current_generate_node.const_pattern());
          new_gar.AddY(literal_info);
          bool has_automorphism = false;
          for (auto& [existed_gar, 
                      existed_gar_lhs_literal_idx] : existed_gars) {
            if (!gar::SameGar(existed_gar, new_gar)) {
              continue;
            }
            has_automorphism = true;
            lhs_disjoint_set.Merge(lhs_literal_idx, 
                       existed_gar_lhs_literal_idx);
            break;
          }
          if (has_automorphism) {
            continue;
          }
          existed_gars.emplace_back(new_gar, lhs_literal_idx);
        }
        existed_gars.clear();

        util::Debug("lhs literal number before remove automorphism: " 
                 + std::to_string(lhs_literal_info_not_in_rhs_set.size()));

        std::vector<size_t> lhs_literal_id_set_to_evaluate;
        for (size_t lhs_literal_idx = 0; 
                    lhs_literal_idx < lhs_literal_info_not_in_rhs_set.size(); 
                    lhs_literal_idx++) {
          size_t group_id = lhs_disjoint_set.Find(lhs_literal_idx);
          lhs_literal_id_set_to_evaluate.emplace_back(group_id);
        }
        std::sort( lhs_literal_id_set_to_evaluate.begin(), 
                   lhs_literal_id_set_to_evaluate.end() );

        lhs_literal_id_set_to_evaluate.erase(
            std::unique( lhs_literal_id_set_to_evaluate.begin(), 
                         lhs_literal_id_set_to_evaluate.end() ), 
                         lhs_literal_id_set_to_evaluate.end() );

        util::Debug("lhs literal number after remove automorphism: " 
                 + std::to_string(lhs_literal_id_set_to_evaluate.size()));

        std::vector<bool> lhs_literal_info_legal(
                          lhs_literal_info_not_in_rhs_set.size(), false);

        for (size_t lhs_literal_idx = 0; 
                    lhs_literal_idx < lhs_literal_id_set_to_evaluate.size(); 
                    lhs_literal_idx++) {
          size_t group_id = lhs_literal_id_set_to_evaluate[lhs_literal_idx] ;
          assert(group_id < lhs_literal_info_not_in_rhs_set.size());
          assert(group_id >= 0);
          const auto& literal_info = lhs_literal_info_not_in_rhs_set[group_id];
          const std::vector<DataGraphLiteralType> literal_set = {literal_info};

          uint64_t supp = 0;
          assert(candidate_set_for_data_graph_set .size()
                            == ctx.data_graph_set_.size());

          assert(candidate_set_removed_nec_for_data_graph_set .size()
                                        == ctx.data_graph_set_.size());

          for (size_t data_graph_idx = 0; 
                      data_graph_idx < ctx.data_graph_set_.size(); 
                      data_graph_idx++){

            auto& data_graph = ctx.data_graph_set_
                               .at(data_graph_idx).data_graph();

            assert(data_graph_idx < candidate_set_for_data_graph_set.size());
            const auto& candidate_set 
                      = candidate_set_for_data_graph_set.at(data_graph_idx);
            if (candidate_set.empty()) {
              assert(candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx).empty());
              continue;
            }

            assert(data_graph_idx < candidate_set_removed_nec_for_data_graph_set.size());
            const auto& candidate_set_removed_nec 
                      = candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx);
            if (candidate_set_removed_nec.empty()) {
              assert(false);
              continue;
            }

            auto [temp_supp, temp_conf] 
               = HasMatch(current_generate_node.pattern(),
                          literal_set,
                          literal_info,
                          data_graph,
                          ctx.normal_to_ml_edge_label_,
                          candidate_set, 
                          candidate_set_removed_nec, 1, true,
                          ctx.graph_basic_statistics_,
                          ctx.restriction_,
                          ctx.time_limit_,
                          ctx.time_limit_per_supp_);

            assert(temp_supp >= 0 && temp_supp <= 1);
            if (temp_supp == 0) {
              continue;
            }
            assert(temp_supp == 1);
            supp = temp_supp;
            break;
          }

          if (supp == 1) {
            assert(!lhs_literal_info_legal[group_id]);
            // not satisfy support bound but is legal
            lhs_literal_info_legal[group_id] = true;
          }
        }
        
        for (size_t lhs_literal_idx = 0; 
                    lhs_literal_idx < lhs_literal_info_not_in_rhs_set.size(); 
                    lhs_literal_idx++) {
          size_t group_id = lhs_disjoint_set.Find(lhs_literal_idx);
          const auto& lhs_literal_info
                    = lhs_literal_info_not_in_rhs_set[lhs_literal_idx];
          if (!lhs_literal_info_legal[group_id]) {
            // is not legal
            continue;
          }
          // is legal
          current_expand_node.AddLhsLiteral(lhs_literal_info_not_in_rhs_set[lhs_literal_idx]);
        }
      }

      util::Debug("current_expand_node.const_rhs_literals().size(): "
                + std::to_string(current_expand_node.const_rhs_literals().size()));
      for (const auto& rhs_literal : current_expand_node.const_rhs_literals()) {
        std::string rhs_literal_str;
        rhs_literal_str << rhs_literal;
        util::Debug("\t" +  rhs_literal_str);
      }

      if (!ctx.restriction_.specified_confidence_bound()) {
        assert(current_expand_node_rhs_literal_conf.empty());
               current_expand_node_rhs_literal_conf.resize(
               current_expand_node_rhs_literal.size(), -1.0);
      }
      assert(current_expand_node_rhs_literal_conf.size()
          == current_expand_node_rhs_literal.size());

      // enumerate the remain literal as the root of literal trees

      if (ctx.store_match_) {
        current_generate_node.InitLiteralTreesWithSuppAndName(current_expand_node_rhs_literal,
                                               current_expand_node_rhs_literal_conf,
                                               current_expand_node_rhs_literal_supp,
                                               current_expand_node_rhs_literal_name,
                                               current_expand_node_rhs_literal_probability);
      } else {
        current_generate_node.InitLiteralTreesWithSupp(current_expand_node_rhs_literal,
                                               current_expand_node_rhs_literal_conf,
                                               current_expand_node_rhs_literal_supp);
      }

      for (size_t literal_idx = 0; 
                  literal_idx < current_expand_node_rhs_literal.size(); 
                  literal_idx++) {
        const auto& rhs_literal = current_expand_node_rhs_literal[literal_idx];
        if (ctx.restriction_.specified_confidence_bound()) {
          if (current_expand_node_rhs_literal_conf[literal_idx] 
           >= ctx.restriction_.confidence_bound()) {
            continue;
          }
        }
        current_expand_node.AddRhsLiteral(rhs_literal);
      }

      for (const auto& rhs_literal
                     : rhs_literal_info_satisfy_supp_bound_set_duplicated) {
        current_expand_node.AddRhsLiteral(rhs_literal);
      }

      if (current_generate_node.LiteralTreesBegin()
       == current_generate_node.LiteralTreesEnd()) {
        // literal trees are empty
        continue;
      }
      
      // expand the literal trees for given layers
      for (auto literal_tree_it  = current_generate_node.LiteralTreesBegin();
                literal_tree_it != current_generate_node.LiteralTreesEnd();
                literal_tree_it++) {
        //std::cout << "new generation tree" << std::endl;
        auto& literal_tree = *literal_tree_it;
        // different from the paper, it is slightly optimized here
        // where each path from root to vertex denotes all literal
        // that is contained in that vertes
        std::vector<DataGraphLiteralType> existed_literals;

        /// <current_expanding_literal_vertex_ptr, boolean flag>
        //   flag, true for expanding and false for recycling
        std::stack<std::pair<LiteralTreePtr, bool>> literal_stack;
        std::set<LiteralTreePtr> visited;

        LiteralTreePtr literal_tree_root 
                     = literal_tree.VertexBegin();

        assert(literal_tree.CountVertex() == 1);
        visited.emplace(literal_tree_root);
        literal_stack.emplace(literal_tree_root, true);

        if (ctx.restriction_.specified_confidence_bound()) {
          const float kConfidence
              = literal_tree_root
              ->template const_attribute<float>(kConfKey);
          if (kConfidence >= ctx.restriction_.confidence_bound()) {
            // satisfied confidence bound, does not 
            // need to consider
            continue;
          }
        }

        assert(literal_tree_root->FindAttribute(kLiteralKey));
        DataGraphLiteralType y_literal 
          = literal_tree_root
          ->template const_attribute<DataGraphLiteralType>(kLiteralKey);

        typename GUNDAM::VertexID<LiteralTreeType>::type vertex_id_allocator = 1;
        typename GUNDAM::  EdgeID<LiteralTreeType>::type   edge_id_allocator = 0;

        if (!ctx.match_gar_using_set_) {

          //std::cout << "match one by one" << std::endl;
          while (!literal_stack.empty()) {
            auto [ current_literal_handle,
                  current_literal_is_first_considered ]
                        = literal_stack.top();
            literal_stack.pop();
          
            // should not be null
            assert(current_literal_handle);
            // label represents level, starts at 1
            assert(current_literal_handle->label() >= 1);
            if (!current_literal_is_first_considered) {
              // is not first considered, recycle
              // the label of the vertex in the literal tree represents
              // the level of it
              if (existed_literals.size() != current_literal_handle->label()) {
                assert(false);
                util::Error("# existed_literals.size() != current_literal_handle->label() !!!");
                util::Error("# existed_literals.size(): " + std::to_string(existed_literals.size())
                          + "# current_literal_handle->label(): " + std::to_string(current_literal_handle->label()));
              }
              assert(existed_literals.size() == current_literal_handle->label());
              existed_literals.pop_back();
              continue;
            }
            // add back to stack, mark as is not considered at the first time
            literal_stack.emplace(current_literal_handle, false);

            assert(current_literal_handle->FindAttribute(kLiteralKey));
            assert(current_literal_handle->FindAttribute(kConfKey));
            // expanding
            existed_literals.emplace_back(
                current_literal_handle
                      ->template const_attribute<DataGraphLiteralType>(kLiteralKey));

            assert(!ctx.restriction_.specified_confidence_bound()
                || (current_literal_handle
                  ->template const_attribute<float>(kConfKey) < ctx.restriction_.confidence_bound()) );
                        
            if (current_literal_handle->label() == ctx.j_) {
              assert(ctx.j_ == 1); // otherwise, it should not be added into the stack
              // reached the literal number limitation
              // does not need to expand
              continue;
            }
            assert(current_literal_handle->label() < ctx.j_);

            // the label of the vertex in the literal tree represents
            // the level of it
            if (existed_literals.size() != current_literal_handle->label()) {
              util::Error("* existed_literals.size() != current_literal_handle->label()!");
              util::Error("* existed_literals.size(): " + std::to_string(existed_literals.size()));
              util::Error("* current_literal_handle->label(): " + std::to_string(current_literal_handle->label()));
              assert(false);
            }
            assert(existed_literals.size() == current_literal_handle->label());

            std::vector<DataGraphLiteralType> expandable_literals;

            // expandable_literals = current_expand_node.lhs_literals
            //                                     - existed_literals
            // legal literals except the ones have been considered
            ExpandableLhsLiterals(current_expand_node.const_lhs_literals(), 
                                                        existed_literals,
                                                    expandable_literals);

            // used to find automorphism gars
            std::vector<GarType> existed_gars;
            for (const auto& expandable_literal : expandable_literals) {
              // edge literal would not be added in lhs
              assert(expandable_literal.literal_type() 
                        != gar::LiteralType::kEdgeLiteral);
              existed_literals.emplace_back(expandable_literal);
              assert(existed_literals.size() > 1);
              if (ctx.restriction_.literals_connected()) {
                std::set<VertexIDType> vertexes_id_in_literals;
                for (const auto& existed_literal : existed_literals) {
                  for (const auto& vertex_id : existed_literal.vertex_id_set()) {
                    vertexes_id_in_literals.emplace(vertex_id);
                  }
                }
                assert(!vertexes_id_in_literals.empty());
                if (!GUNDAM::Connected(current_generate_node.const_pattern(),
                                      vertexes_id_in_literals)) {
                  existed_literals.pop_back();
                  continue;
                }
              }

              GarType new_gar(current_generate_node.const_pattern());
              for (const auto& existed_literal : existed_literals) {
                new_gar.AddX(existed_literal);
              }
              bool has_same = false;
              for (const auto& existed_gar : existed_gars){
                if (gar::SameGar(existed_gar, new_gar)){
                  has_same = true;
                  break;
                }
              }
              if (has_same) {
                existed_literals.pop_back();
                continue;
              }
              existed_gars.emplace_back(std::move(new_gar));

              assert(candidate_set_for_data_graph_set .size()
                                == ctx.data_graph_set_.size());

              assert(candidate_set_removed_nec_for_data_graph_set .size()
                                            == ctx.data_graph_set_.size());

              // whether has satisfied support bound on at least
              // one data graph
              bool satisfy_supp_bound = false,
                  satisfy_conf_bound = false;

              float conf = -1.0;
              uint64_t supp = 0;
              double probability = 0.0;
              std::vector<std::pair<std::vector<VertexIDType>,
                                    std::vector<EdgeIDType>>> match_vec;

              for (size_t data_graph_idx = 0; 
                          data_graph_idx < ctx.data_graph_set_.size(); 
                          data_graph_idx++) {
                auto& data_graph = ctx.data_graph_set_.at(data_graph_idx).data_graph();

                assert(data_graph_idx < candidate_set_for_data_graph_set.size());
                const auto& candidate_set 
                          = candidate_set_for_data_graph_set.at(data_graph_idx);
                if (candidate_set.empty()) {
                  assert(candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx).empty());
                  continue;
                }

                assert(data_graph_idx < candidate_set_removed_nec_for_data_graph_set.size());
                const auto& candidate_set_removed_nec 
                          = candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx);
                if (candidate_set_removed_nec.empty()) {
                  assert(false);
                  continue;
                }

                auto [temp_supp, temp_conf]
                  = HasMatch(current_generate_node.pattern(),
                              existed_literals,
                            *existed_literals.begin(),
                              data_graph,
                              ctx.normal_to_ml_edge_label_,
                              candidate_set,
                              candidate_set_removed_nec,
                              ctx.support_bound_, false,
                              ctx.graph_basic_statistics_,
                              ctx.restriction_,
                              ctx.time_limit_,
                              ctx.time_limit_per_supp_,
                              ctx.store_match_,
                              match_vec,
                              probability);
                //std::cout << "line 2849 prob " << probability << std::endl;
                
                if (temp_supp > supp) {
                  supp = temp_supp;
                }
                if (temp_supp < ctx.support_bound_) {
                  // does not satisfy support bound
                  continue;
                }
                satisfy_supp_bound = true;
                if (!ctx.restriction_.specified_confidence_bound()) {
                  break;
                }
                if (temp_conf >= ctx.restriction_.confidence_bound()) {
                  satisfy_conf_bound = true;
                  conf = temp_conf;
                  break;
                }
              }

              // this literal has already been considered, maintain
              // the existed_literals for other literals in this level
              existed_literals.pop_back();

              // for all expandable literal, add a vertex in the literal
              // tree that contains it
              auto [ add_vertex_handle,
                    add_vertex_ret ] = literal_tree.AddVertex(vertex_id_allocator++, 
                                                              current_literal_handle->label() + 1);
              // added successfully
              assert(add_vertex_handle);
              assert(add_vertex_ret);
              // label needs to be legal
              assert(add_vertex_handle->label() <= ctx.j_);

              auto [ add_edge_handle,
                    add_edge_ret ] = literal_tree.AddEdge(current_literal_handle->id(), 
                                                                add_vertex_handle->id(),
                                                          kLiteralTreeDefaultEdgeLabel, 
                                                          edge_id_allocator++);
              // added successfully
              assert(add_edge_handle);
              assert(add_edge_ret);

              auto [ attr_handle,
                    attr_ret ] = add_vertex_handle->AddAttribute(kLiteralKey, 
                                                                  expandable_literal);
              // added successfully
              assert(attr_ret);

              auto [ supp_attr_handle,
                      supp_attr_ret ] = add_vertex_handle->AddAttribute(kSuppKey, supp);
              assert(supp_attr_handle);
              assert(supp_attr_ret);
              //added successfully


              auto [ conf_attr_handle,
                    conf_attr_ret ] = add_vertex_handle->AddAttribute(kConfKey,conf);
              // added successfully
              assert(conf_attr_ret);


              if (ctx.store_match_) {
                std::string gar_name = ExportGARMatch(*existed_literals.begin(),
                                                  match_vec, kCurrentPatternEdgeSize,
                                                  process_expand_node_idx, gar_offset,
                                                  ctx.fid_, ctx.output_match_dir_);
                gar_offset++;
                auto [name_attr_handle,
                      name_attr_ret ] = add_vertex_handle->AddAttribute(kMatchFileName, gar_name);
                auto [prob_attr_handle,
                      prob_attr_ret] = add_vertex_handle->AddAttribute(kProbability, probability);
              }

              if (!satisfy_supp_bound) {
                // this gar does not satisfy support bound on any data
                // graph, does not needs to be considered
                continue;
              }

              if (add_vertex_handle->label() == ctx.j_) {
                // has reached the literal number limitation,
                // does not need to add int the stack
                continue;
              }
              if (ctx.restriction_.specified_confidence_bound()
              && satisfy_conf_bound) {
                // does no need to be further expanded
                continue;
              }
              literal_stack.emplace(add_vertex_handle, true);
            }
          }
        } else {
          //std::cout << "match by set" << std::endl;
          while (!literal_stack.empty()) {
            //std::cout << "new literal in literal stack to check" << std::endl;
            auto [ current_literal_handle,
                  current_literal_is_first_considered ]
                        = literal_stack.top();
            literal_stack.pop();
          
            // should not be null
            assert(current_literal_handle);
            // label represents level, starts at 1
            assert(current_literal_handle->label() >= 1);
            if (!current_literal_is_first_considered) {
    
            //std::cout << "here after check0.5" << std::endl;
              // is not first considered, recycle
              // the label of the vertex in the literal tree represents
              // the level of it
              if (existed_literals.size() != current_literal_handle->label()) {
                std::cout << "error" << std::endl;
                assert(false);
                util::Error("# existed_literals.size() != current_literal_handle->label() !!!");
                util::Error("# existed_literals.size(): " + std::to_string(existed_literals.size())
                          + "# current_literal_handle->label(): " + std::to_string(current_literal_handle->label()));
              }
              assert(existed_literals.size() == current_literal_handle->label());
              //std::cout << "existed_literals size" << existed_literals.size() << std::endl;
              existed_literals.pop_back();
              //std::cout << "existed_literals size after" << existed_literals.size() << std::endl;
              //std::cout << "stack size " << literal_stack.size() << std::endl;
              continue;
            }
            //std::cout << "here after check" << std::endl;
            // add back to stack, mark as is not considered at the first time
            literal_stack.emplace(current_literal_handle, false);

            assert(current_literal_handle->FindAttribute(kLiteralKey));
            assert(current_literal_handle->FindAttribute(kConfKey));
            //std::cout << "here after check 1" << std::endl;
            // expanding
            existed_literals.emplace_back(
                current_literal_handle
                      ->template const_attribute<DataGraphLiteralType>(kLiteralKey));

            assert(!ctx.restriction_.specified_confidence_bound()
                || (current_literal_handle
                  ->template const_attribute<float>(kConfKey) < ctx.restriction_.confidence_bound()) );
                        
            if (current_literal_handle->label() == ctx.j_) {
              assert(ctx.j_ == 1); // otherwise, it should not be added into the stack
              // reached the literal number limitation
              // does not need to expand
              continue;
            }
            assert(current_literal_handle->label() < ctx.j_);
            //std::cout << "here after check 2" << std::endl;

            // the label of the vertex in the literal tree represents
            // the level of it
            if (existed_literals.size() != current_literal_handle->label()) {
              util::Error("* existed_literals.size() != current_literal_handle->label()!");
              util::Error("* existed_literals.size(): " + std::to_string(existed_literals.size()));
              util::Error("* current_literal_handle->label(): " + std::to_string(current_literal_handle->label()));
              assert(false);
            }
            assert(existed_literals.size() == current_literal_handle->label());

            std::vector<DataGraphLiteralType> expandable_literals;
            //std::cout << "begin checking expandbale literals" << std::endl;
            // expandable_literals = current_expand_node.lhs_literals
            //                                     - existed_literals
            // legal literals except the ones have been considered
            ExpandableLhsLiterals(current_expand_node.const_lhs_literals(), 
                                                        existed_literals,
                                                    expandable_literals);

            // used to find automorphism gars
            std::vector<GarType> existed_gars;

            //std::cout << "expandale lhs literal size " << expandable_literals.size() << std::endl;
            
            std::vector<DataGraphLiteralType> new_literals_to_check;
            GarType common_gar(current_generate_node.const_pattern());
            for (const auto& existed_literal : existed_literals) {
              common_gar.AddX(existed_literal);
            }

            //std::cout << "here expandable start" << std::endl;
            for (const auto& expandable_literal : expandable_literals) {
              //std::cout << "check" << std::endl;
              // edge literal would not be added in lhs
              assert(expandable_literal.literal_type() 
                        != gar::LiteralType::kEdgeLiteral);
              existed_literals.emplace_back(expandable_literal);
              assert(existed_literals.size() > 1);
              if (ctx.restriction_.literals_connected()) {
                std::set<VertexIDType> vertexes_id_in_literals;
                for (const auto& existed_literal : existed_literals) {
                  for (const auto& vertex_id : existed_literal.vertex_id_set()) {
                    vertexes_id_in_literals.emplace(vertex_id);
                  }
                }
                assert(!vertexes_id_in_literals.empty());
                if (!GUNDAM::Connected(current_generate_node.const_pattern(),
                                      vertexes_id_in_literals)) {
                  existed_literals.pop_back();
                  continue;
                }
              }

              GarType new_gar(current_generate_node.const_pattern());
              for (const auto& existed_literal : existed_literals) {
                new_gar.AddX(existed_literal);
              }
              bool has_same = false;
              for (const auto& existed_gar : existed_gars){
                if (gar::SameGar(existed_gar, new_gar)){
                  has_same = true;
                  break;
                }
              }
              if (has_same) {
                existed_literals.pop_back();
                continue;
              }
              existed_gars.emplace_back(std::move(new_gar));

              assert(candidate_set_for_data_graph_set .size()
                                == ctx.data_graph_set_.size());

              assert(candidate_set_removed_nec_for_data_graph_set .size()
                                            == ctx.data_graph_set_.size());

              new_literals_to_check.emplace_back(expandable_literal);


              // this literal has already been considered, maintain
              // the existed_literals for other literals in this level
              existed_literals.pop_back();
            }
            //std::cout << "new _literals size" << new_literals_to_check.size() << std::endl;
            if (new_literals_to_check.size() == 0) continue;
            
            std::vector<bool> satisfy_supp_bound(new_literals_to_check.size(), false),
                              satisfy_conf_bound(new_literals_to_check.size(), false);
            // // whether has satisfied support bound on at least
            // // one data graph
            std::vector<float> conf(new_literals_to_check.size(), -1.0);
            std::vector<uint64_t> supp(new_literals_to_check.size(), 0);


            // bool satisfy_supp_bound = false,
            //       satisfy_conf_bound = false;

            // float conf = -1.0;
            // int supp = 0;
            std::vector<std::vector<std::pair<std::vector<VertexIDType>, std::vector<EdgeIDType>>>>
                match_vec_for_all;
            for (size_t data_graph_idx = 0; 
                        data_graph_idx < ctx.data_graph_set_.size(); 
                        data_graph_idx++) {
              auto& data_graph = ctx.data_graph_set_.at(data_graph_idx).data_graph();

              assert(data_graph_idx < candidate_set_for_data_graph_set.size());
              const auto& candidate_set 
                        = candidate_set_for_data_graph_set.at(data_graph_idx);
              if (candidate_set.empty()) {
                assert(candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx).empty());
                continue;
              }

              assert(data_graph_idx < candidate_set_removed_nec_for_data_graph_set.size());
              const auto& candidate_set_removed_nec 
                        = candidate_set_removed_nec_for_data_graph_set.at(data_graph_idx);
              if (candidate_set_removed_nec.empty()) {
                assert(false);
                continue;
              }


              //std::cout << "start to match set" << std::endl;
              auto [temp_supp_vec, temp_conf_vec]
                = HasMatchForSet(current_generate_node.pattern(),
                            existed_literals,
                            *existed_literals.begin(),
                            new_literals_to_check,
                            data_graph,
                            ctx.normal_to_ml_edge_label_,
                            candidate_set,
                            candidate_set_removed_nec,
                            ctx.support_bound_, false,
                            ctx.graph_basic_statistics_,
                            ctx.restriction_,
                            ctx.time_limit_,
                            ctx.time_limit_per_supp_,
                            ctx.store_match_,
                            match_vec_for_all);

              if (temp_supp_vec.size() == 0) {
                std::cout << "error occurred" << std::endl;
              }

             //for (unsigned i = 0; i < temp_supp_vec.size(); i++) {
             //   std::cout << "support " << temp_supp_vec[i] << std::endl;
             //}

              for (unsigned new_literal_idx = 0;
                            new_literal_idx < new_literals_to_check.size();
                            new_literal_idx++) {
                if (temp_supp_vec[new_literal_idx] >= supp[new_literal_idx]) {
                  supp[new_literal_idx] = temp_supp_vec[new_literal_idx];
                }
                if (temp_supp_vec[new_literal_idx] < ctx.support_bound_) {
                  // does not satisfy support bound
                  satisfy_supp_bound[new_literal_idx] = false;
                } else {
                  satisfy_supp_bound[new_literal_idx] = true;
                }

                if (!ctx.restriction_.specified_confidence_bound()) {
                  continue;
                }
                if (temp_conf_vec[new_literal_idx] >= ctx.restriction_.confidence_bound()) {
                  satisfy_conf_bound[new_literal_idx] = true;
                  conf[new_literal_idx] = temp_conf_vec[new_literal_idx];
                  //break;
                }
              }
            }
            //std::cout << "here 4" << std::endl;

            for (unsigned new_literal_idx = 0;
                            new_literal_idx < new_literals_to_check.size();
                            new_literal_idx++) {

              // for all expandable literal, add a vertex in the literal
              // tree that contains it
              auto [ add_vertex_handle,
                      add_vertex_ret ] = literal_tree.AddVertex(vertex_id_allocator++, 
                                                                current_literal_handle->label() + 1);
              // added successfully
              assert(add_vertex_handle);
              assert(add_vertex_ret);
              // label needs to be legal
              assert(add_vertex_handle->label() <= ctx.j_);

              if (!add_vertex_ret) {
                std::cout << "cant add vertex " << std::endl;
              }

              auto [ add_edge_handle,
                      add_edge_ret ] = literal_tree.AddEdge(current_literal_handle->id(), 
                                                                add_vertex_handle->id(),
                                                          kLiteralTreeDefaultEdgeLabel, 
                                                            edge_id_allocator++);
              // added successfully
              assert(add_edge_handle);
              assert(add_edge_ret);

                            if (!add_edge_ret) {
                std::cout << "cant add edge " << std::endl;
              }

              auto [ attr_handle,
                      attr_ret ] = add_vertex_handle->AddAttribute(kLiteralKey, 
                                                            new_literals_to_check[new_literal_idx]);
              // added successfully
              assert(attr_ret);
              if (!attr_ret) {
                std::cout << "cant add attr " << std::endl;
              }

              auto [ conf_attr_handle,
                      conf_attr_ret ] = add_vertex_handle->AddAttribute(kConfKey, conf[new_literal_idx]);
              // added successfully
              assert(conf_attr_ret);

              if (!conf_attr_ret) {
                std::cout << "cant add conf " << std::endl;
              }

              //if (supp[new_literal_idx] == 0) {
              //  std::cout << "support is 0" << std::endl;
              //}
              //auto [ supp_attr_handle,
              //        supp_attr_ret ] = add_vertex_handle->AddAttribute(kSuppKey, supp[new_literal_idx]);
              auto [ supp_attr_handle,
                      supp_attr_ret ] = add_vertex_handle->AddAttribute(kSuppKey, supp[new_literal_idx]);
              assert(supp_attr_handle);
              assert(supp_attr_ret);
              //added successfully
                            if (!supp_attr_ret) {
                std::cout << "cant add support " << std::endl;
              }


              if (ctx.store_match_) {
                std::string gar_name = ExportGARMatch(*existed_literals.begin(),
                                                      match_vec_for_all[new_literal_idx],
                                                      kCurrentPatternEdgeSize,
                                                      process_expand_node_idx,
                                                      gar_offset, ctx.fid_, 
                                                      ctx.output_match_dir_);
                gar_offset++;
                auto [name_attr_handle,
                      name_attr_ret ] = add_vertex_handle->AddAttribute(kMatchFileName, gar_name);
              }

              if (!satisfy_supp_bound[new_literal_idx]) {
                // this gar does not satisfy support bound on any data
                // graph, does not needs to be considered
                continue;
              }

              if (add_vertex_handle->label() == ctx.j_) {
                // has reached the literal number limitation,
                // does not need to add int the stack
                continue;
              }
              if (ctx.restriction_.specified_confidence_bound()
                && satisfy_conf_bound[new_literal_idx]) {
                // does no need to be further expanded
                continue;
              }
              literal_stack.emplace(add_vertex_handle, true);
            }
            //std::cout << "here 5" << std::endl;
          }
        }
      }
      omp_set_lock(&current_generate_level_lock);
      current_generate_level.AddGenerateTreeNode(std::move(current_generate_node));
      omp_unset_lock(&current_generate_level_lock);
    }

    std::cout << "start to export on level " << std::endl;

    ExportGeneratedGarsWithSuppAndTree<GraphPatternType, 
                           DataGraphType>(
        kCurrentPatternEdgeSize, 
          current_generate_level,
          ctx.restriction_,
          ctx.support_bound_,
          ctx.output_gar_dir_,
          ctx.fid_);

    // ExportGeneratedGars<GraphPatternType, 
    //                        DataGraphType>(
    //     kCurrentPatternEdgeSize, 
    //       current_generate_level,
    //       ctx.restriction_,
    //       ctx.output_gar_dir_,
    //       ctx.fid_);

    auto& channel_0 = messages.Channels()[0];
    // send expand_nodes, where literal were updated, back to kExpandPatternFragID
    std::string msg(kUpdateExpandNodesPrefix);
    if (ctx.fid_ != kExpandPatternFragID) {
      // the kExpandPatternFragID does not need to transform data through the
      // MPI method, 
      assert(ctx.expand_level_  .size() 
          == ctx.process_expand_node_id_set_.size());
      for (auto& expand_node : ctx.expand_level_) {
        UpdateLiteral<DataGraphLiteralType> update_literal;
        update_literal.set_id(expand_node.id());
        update_literal.rhs_literals().swap(expand_node.rhs_literals());
        update_literal.lhs_literals().swap(expand_node.lhs_literals());
        msg << update_literal;
      }
      ctx.expand_level_.clear();
    }
    channel_0.SendToFragment(kExpandPatternFragID, msg);

#ifdef PROFILING
    ctx.preprocess_time -= GetCurrentTime();
#endif

#ifdef PROFILING
    ctx.preprocess_time += GetCurrentTime();
    ctx.exec_time -= GetCurrentTime();
#endif

#ifdef PROFILING
    ctx.exec_time += GetCurrentTime();
    ctx.postprocess_time -= GetCurrentTime();
#endif

#ifdef PROFILING
    ctx.postprocess_time += GetCurrentTime();
#endif

    return;
  }
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GAR_DISCOVER_H_
