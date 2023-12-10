#ifndef EXAMPLES_ANALYTICAL_APPS_INC_GAR_DISCOVER_GAR_DISCOVER_H_
#define EXAMPLES_ANALYTICAL_APPS_INC_GAR_DISCOVER_GAR_DISCOVER_H_

#include <grape/grape.h>
#include <omp.h>
#include <chrono>

#include "../fragment2gundam.h"
#include "../fragment_graph.h"
#include "../fragment_graph_with_index.h"
#include "../timer.h"

#include "util/load_rule.h"
#include "rule_discover/inc_gar_discover_context.h"
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
class IncGarDiscover : public ParallelAppBase<FRAG_T, IncGarDiscoverContext<FRAG_T>>,
                    public ParallelEngine {
 public:
  static constexpr LoadStrategy load_strategy = LoadStrategy::kLoadWholeGraph;

 private:
  using     VertexIDType = typename IncGarDiscoverContext<FRAG_T>::    VertexIDType;
  using  VertexLabelType = typename IncGarDiscoverContext<FRAG_T>:: VertexLabelType;
  using       EdgeIDType = typename IncGarDiscoverContext<FRAG_T>::      EdgeIDType;
  using    EdgeLabelType = typename IncGarDiscoverContext<FRAG_T>::   EdgeLabelType;
  using    DataGraphType = typename IncGarDiscoverContext<FRAG_T>::   DataGraphType;
  using GraphPatternType = typename IncGarDiscoverContext<FRAG_T>::GraphPatternType;
  
  using DataGraphLiteralType 
          = gar::LiteralInfo<GraphPatternType,
                                DataGraphType>;

  using VertexAttributeKeyType = typename IncGarDiscoverContext<FRAG_T>::VertexAttributeKeyType;

  using LiteralStandAloneInfoType 
 = gar::LiteralStandAloneInfo<GraphPatternType,
                                 DataGraphType>;

  using GarType = gar::GraphAssociationRule<GraphPatternType,
                                               DataGraphType>;
  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using IncGraphPackage = _gar_discover::IncGraphPackage<DataGraphType>;
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
  INSTALL_PARALLEL_WORKER(IncGarDiscover<FRAG_T>, 
                          IncGarDiscoverContext<FRAG_T>,
                          FRAG_T)
  using vertex_t = typename fragment_t::vertex_t;

  /**
   * @brief Partial evaluation for IncGarDiscover.
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

    using GenerateTreeNodeType 
        = GenerateTreeNode<GraphPatternType,
                              DataGraphType>;


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

    std::tuple<std::string,
                           std::string,
                           std::string> data_graph_path;

    auto [ data_graph_string, 
            data_graph_ret ] = util::GetDataGraphInfoFromYaml(data_graph_path_config);
    if (!data_graph_ret) {
      util::Error("load data graph error!");
      return;
    }
    data_graph_path = data_graph_string;


    if (!config["RulePath"]) {
      util::Error("cannot get rule path");
      return;
    }
    YAML::Node rule_path_config = config["RulePath"];

    //std::vector<std::pair<GarType, std::string>> rule_set;
    std::unordered_map<int, int> rule_id2parent_id;
    std::unordered_map<int, typename DataGraphType::VertexCounterType> rule_id2supp;
    std::unordered_map<int, std::string> rule_id2pivot_file;
    std::unordered_map<int, double> rule_id2probability;
    if (rule_path_config.IsSequence()) {
      util::Error("should be a gar set file");
      return;
    }
    else {
      ctx.input_gar_set_ = util::LoadRuleWithTreeAndSupp<
                        GraphPatternType, DataGraphType>(rule_path_config,
                                            rule_id2parent_id, rule_id2supp,
                                            rule_id2pivot_file, rule_id2probability);
      std::cout << "gar set size " << ctx.input_gar_set_.size() << std::endl;
    }

    ctx.rule_id2parent_id_ = rule_id2parent_id;
    ctx.rule_id2supp_ = rule_id2supp;


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
    // ##   load data graph  through the loaded path           ##
    // ##########################################################

    const auto& data_graph_v_file = std::get<0>(data_graph_path);
    const auto& data_graph_e_file = std::get<1>(data_graph_path);
    const auto& ml_literal_edge_file = std::get<2>(data_graph_path);

    ctx.data_graph_.emplace_back(specified_edge_type_set,
                                      specified_edge_label_set,
                                      ctx.restriction_.constant_freq_bound(),
                                      data_graph_v_file,
                                      data_graph_e_file,
                                      ctx.origin_graph_basic_statistics_,
                                      ctx.updated_graph_basic_statistics_,
                                      ctx.knowledge_graph_v_file_,
                                      ctx.knowledge_graph_e_file_,
                                      ctx.er_file_);
    if (!ctx.data_graph_[0].load_data_graph_success()) {
      util::Error("load inc data graph fail!");
      return;
    }
    ctx.data_graph_[0].GenerateGraphNec();

    auto& origin_data_graph = ctx.data_graph_[0].origin_data_graph();
    auto& updated_data_graph = ctx.data_graph_[0].updated_data_graph();
    auto& origin_data_graph_nec = ctx.data_graph_[0].origin_data_graph_nec();
    auto& updated_data_graph_nec = ctx.data_graph_[0].updated_data_graph_nec();
    auto& id_to_origin_updated_vertex_handle = ctx.data_graph_[0].id_to_vertex_handle();
    auto& updated_vertex_id = ctx.data_graph_[0].updated_vertex_id();

    // #####################################################
    // ##  generate map from edge_label to ml_edge_label  ##
    // ##  as well as ml_edge_label to edge_label         ##
    // #####################################################
    EdgeLabelType max_edge_label = 0,
      minimal_negtive_edge_label = 0;
    for (const auto& [edge_label,
                      edge_label_counter]
                : ctx.origin_graph_basic_statistics_
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
                            : ctx.origin_graph_basic_statistics_
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
                           : ctx.origin_graph_basic_statistics_
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
    const std::string& ml_e_file = std::get<2>(data_graph_path);
  
    if (ml_e_file == "") {
      // does not have ml edge for this data graph
    } else {
      if (!ctx.data_graph_[0].AddMlEdge(ml_e_file,
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

    ctx.threshold_ = 1.0;
    if (config["SupportEstimation"]) {
      ctx.support_estimation_ = config["SupportEstimation"].as<bool>();
    } else {
      ctx.support_estimation_ = false;
    }

    if (config["MatchUsingSet"]) {
      ctx.match_gar_using_set_ = config["MatchUsingSet"].as<bool>();
    } else {
      ctx.match_gar_using_set_ = false;
    }
    

    if (ctx.support_estimation_) {
      ctx.threshold_ = 0.99;
      if (config["EstimationThreshold"]) {
        ctx.threshold_ = config["EstimationThreshold"].as<double>();
      }
    }

    int current_worker_id, current_worker_num;
    MPI_Comm_rank(MPI_COMM_WORLD, &current_worker_id);
    MPI_Comm_size(MPI_COMM_WORLD, &current_worker_num);

    std::cout << "worker id " << current_worker_id << " worker_num " << current_worker_num << std::endl;

    if (!ctx.RebuildLocalGARGenerateTree(
                                      rule_id2parent_id, rule_id2supp,
                                      rule_id2pivot_file, rule_id2probability,
                                      current_worker_id, current_worker_num)) {
      std::cout << "failed to rebuild generate tree!" << std::endl;
      return;
    }

    util::Info("rule loaded, total number: " + std::to_string(ctx.input_gar_set_.size()));


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

    if (config["WithIndex"]) {
      ctx.with_index_ = config["WithIndex"].as<bool>();
      if (ctx.with_index_) {
        std::cout << "######inc gar discover will use index######" << std::endl;
      } else {
        std::cout << "######inc gar discover will not use index######" << std::endl;
      }
    } else {
      std::cout << "Default : ###inc gar discover will not use index####" << std::endl;
      ctx.with_index_ = false;
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
        std::cout << "########## should not be here #####" << std::endl;
        return;

        // // ################################
        // // ##   specified support file   ##
        // // ################################
        // const size_t kLiteralVertexNum = stand_alone_literal_info.VertexNum();
        // if (kLiteralVertexNum != 2) {
        //   util::Error("support only literals contains two vertex now!");
        //   return;
        // }
        // if (stand_alone_literal_info.literal_type() != gar::LiteralType::kEdgeLiteral) {
        //   // variable literal has problem in simulation with specified support set
        //   // in HasMatch
        //   util::Error("support only edge literals contains two vertex now!");
        //   return;
        // }

        // const std::string kLiteralSuppFile 
        //    = specified_rhs_literal_config["SuppFile"].as<std::string>();

        // std::ifstream literal_supp_file(kLiteralSuppFile);

        // if (!literal_supp_file.good()) {
        //   util::Error ( "literal_supp_file of the " 
        //               + std::to_string(specified_rhs_literal_set_config_idx) 
        //               + "'th literal: " + kLiteralSuppFile
        //               + " is not good! ");
        //   return;
        // }

        // std::vector<
        // std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>> literal_supp_set;

        // while (literal_supp_file) {
        //   std::string s;
        //   if (!std::getline( literal_supp_file, s ))
        //     break;

        //   std::istringstream ss( s );
        //   std::vector <std::string> record;
        //   std::string buf; 
        //   while (ss >> buf)
        //     record.emplace_back(buf);

        //   std::vector <VertexIDType> support_id_set;
        //   for (const auto& id_str : record) {
        //     support_id_set.emplace_back(std::stoi(id_str));
        //   }
        //   if (kLiteralVertexNum != support_id_set.size()) {
        //     if ((kLiteralVertexNum + 1) != support_id_set.size()
        //      || (support_id_set[kLiteralVertexNum] != 0
        //       && support_id_set[kLiteralVertexNum] != 1)) {
        //       util::Error("support set miss match in line: " 
        //                  + std::to_string(literal_supp_set.size())
        //                  + " of the " 
        //                  + std::to_string(specified_rhs_literal_set_config_idx) 
        //                  + "'th literal");
        //       return;
        //     }
        //   }
        //   if (kLiteralVertexNum + 1 == support_id_set.size()) {
        //     assert(support_id_set[kLiteralVertexNum] == 0
        //         || support_id_set[kLiteralVertexNum] == 1);
        //     support_id_set.pop_back();
        //   }
        //   std::vector <typename GUNDAM::VertexHandle<DataGraphType>::type> support;
        //   support.reserve(support_id_set.size());
        //   // support set only allowed to be specified on one data graph
        //   for (const auto& vertex_id : support_id_set) {
        //     auto vertex_handle = ctx.data_graph_.
        //                             .origin_data_graph().FindVertex(vertex_id);
        //     if (!vertex_handle) {
        //       util::Error("vertex: " + std::to_string(vertex_id) 
        //                 + " cannot be found in origin data graph!");
        //       return;
        //     }
        //     support.emplace_back(vertex_handle);
        //   }
        //   assert(support.size() == support_id_set.size());
        //   literal_supp_set.emplace_back(support);
        // }
        // if (literal_supp_set.empty()) {
        //   util::Error("Specified empty support set for the "
        //              + std::to_string(specified_rhs_literal_set_config_idx) 
        //              + "'th literal");
        //   return;
        // }
        // std::sort(literal_supp_set.begin(),
        //           literal_supp_set.end());
        // ctx.restriction_.AddSpecifiedRhsLiteral(stand_alone_literal_info,
        //                                         literal_supp_set);
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
    } else {
      ctx.restriction_.SpecifyRuleType("gar");
    }
    assert(ctx.restriction_.gar()
        || ctx.restriction_.gcr()
        || ctx.restriction_.horn_rule());

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

    unsigned pass_count = 0;
    ctx.begin_time_ = std::chrono::high_resolution_clock::now();
    std::vector<int> gar_to_extend;
    omp_lock_t gar_id_lock;
    omp_init_lock(&gar_id_lock);

    omp_set_nested(false);

    #pragma omp parallel for schedule(dynamic) 
    for (size_t generate_node_idx = 0; 
                generate_node_idx < static_cast<size_t>(ctx.generate_tree_level_.size()); 
                //generate_node_idx < 100;
                generate_node_idx++) {
      auto &current_generate_tree_node = ctx.generate_tree_level_.node(generate_node_idx);
      auto &pattern = current_generate_tree_node.pattern();

      std::cout << "pattern generate_node_idx " << generate_node_idx << std::endl;

      if (!ctx.match_gar_using_set_) {
        for (auto literal_tree_it = current_generate_tree_node.LiteralTreesBegin();
                  literal_tree_it != current_generate_tree_node.LiteralTreesEnd();
                  literal_tree_it++) {
          auto &literal_tree = *literal_tree_it;
          typename GUNDAM::VertexID<LiteralTreeType>::type kLiteralTreeRootId = 0;

          //std::cout << "new literal tree " << std::endl;

          using LiteralVertexHandleType
                    = typename GUNDAM::VertexHandle<LiteralTreeType>::type;
          
          auto root_handle = literal_tree.FindVertex(kLiteralTreeRootId);
          assert(root_handle->CountInVertex() == 0);

          std::queue<LiteralVertexHandleType> literal_queue;
          std::map<LiteralVertexHandleType, std::vector<DataGraphLiteralType>> literal_tree_handle2literals;
          DataGraphLiteralType y_literal 
              = root_handle->template const_attribute<DataGraphLiteralType>(_gar_discover::kLiteralKey);

          std::vector<DataGraphLiteralType> root_literals;
          root_literals.emplace_back(y_literal);

          literal_queue.emplace(root_handle);
          literal_tree_handle2literals[root_handle] = root_literals;

          std::map<LiteralVertexHandleType, unsigned> child_num_with_valid_supp;

          while (!literal_queue.empty()) {
            auto current_literal_handle = literal_queue.front();
            literal_queue.pop();

            auto current_literal_new_supp_handle
                    = current_literal_handle->FindAttribute(kNewSuppKey);

            if (current_literal_new_supp_handle) {
              auto current_literal_supp = current_literal_handle
                                ->template const_attribute<int>(_gar_discover::kNewSuppKey);
              if (current_literal_supp < ctx.support_bound_) {
                continue;
              }
            } else {
              //std::cout << "exist no newkey node " << std::endl;
            }

            unsigned valid_child_num = 0;

            for (auto out_edge_it = current_literal_handle->OutEdgeBegin();
                    !out_edge_it.IsDone();
                      out_edge_it++) {
              auto dst_literal_handle = out_edge_it->dst_handle();
              auto dst_supp = dst_literal_handle
                                ->template const_attribute<int>(_gar_discover::kSuppKey);
              if (dst_literal_handle->FindAttribute(kNewSuppKey)) {
                dst_supp = dst_literal_handle
                                ->template const_attribute<int>(_gar_discover::kNewSuppKey);
              }
              if(dst_supp >= ctx.support_bound_) {
                valid_child_num++;
              }
            }
            child_num_with_valid_supp[current_literal_handle] = valid_child_num;


            auto &current_literals = literal_tree_handle2literals[current_literal_handle];
            if (valid_child_num == 0 && (!current_literal_new_supp_handle)) {
              // check the gar with no valid child and not checked yet
              double probability = current_literal_handle
                                ->template const_attribute<double>(_gar_discover::kProbability);
              using CandidateSetContainer 
                    = std::map<GraphPatternVertexHandleType,
                            std::vector<DataGraphVertexHandleType>>;


              CandidateSetContainer origin_candidate_set,
                                  updated_candidate_set,
                                    origin_candidate_set_removed_nec,
                                  updated_candidate_set_removed_nec;

              bool init_origin_candidate_set_succ = false;
              bool init_updated_candidate_set_succ = false;

              init_origin_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                    pattern, origin_data_graph, origin_candidate_set);

              
              init_updated_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                    pattern, updated_data_graph, updated_candidate_set);
              

              if (!init_origin_candidate_set_succ) {
                origin_candidate_set.clear();
              }

              if (!init_updated_candidate_set_succ) {
                updated_candidate_set.clear();
              }

              bool refine_origin_candidate_set_succ = false;
              bool refine_updated_candidate_set_succ = false;


              refine_origin_candidate_set_succ
                = GUNDAM::_dp_iso_using_match
                        ::RefineCandidateSet(
                  pattern, origin_data_graph, origin_candidate_set);


              refine_updated_candidate_set_succ
                = GUNDAM::_dp_iso_using_match
                        ::RefineCandidateSet(
                  pattern, updated_data_graph, updated_candidate_set);


              if (!refine_origin_candidate_set_succ) {
                origin_candidate_set.clear();
              }

              if (!refine_updated_candidate_set_succ) {
                updated_candidate_set.clear();
              }

              origin_candidate_set_removed_nec = origin_candidate_set;
              GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                  DataGraphType>(
                        origin_candidate_set_removed_nec,
                        origin_data_graph_nec);

              updated_candidate_set_removed_nec = updated_candidate_set;
              GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                  DataGraphType>(
                        updated_candidate_set_removed_nec,
                        updated_data_graph_nec);

              std::string pivot_file_name;
              
              if (ctx.with_index_) {
                auto pivot_file_handle = current_literal_handle
                          ->FindAttribute(_gar_discover::kMatchFileName);
                if (!pivot_file_handle) {
                  std::cout << "can not find pivot file name " << std::endl;
                }
                pivot_file_name = current_literal_handle
                          ->template const_attribute<std::string>(_gar_discover::kMatchFileName);
              }
              auto origin_supp = current_literal_handle
                                  ->template const_attribute<int64_t>(_gar_discover::kSuppKey);
              int64_t estimated_support = origin_supp;
              bool no_need_to_check = IncHasMatchByEstimate(pattern,
                                            current_literals,
                                            root_literals[0],
                                            origin_data_graph,
                                            origin_data_graph_nec,
                                            updated_data_graph,
                                            updated_data_graph_nec,
                                            ctx.normal_to_ml_edge_label_,
                                            origin_candidate_set,
                                            origin_candidate_set_removed_nec,
                                            updated_candidate_set,
                                            updated_candidate_set_removed_nec,
                                            id_to_origin_updated_vertex_handle,
                                            updated_vertex_id,
                                            pivot_file_name,
                                            origin_supp,
                                            ctx.support_bound_,
                                            probability,
                                            ctx.threshold_,
                                            estimated_support,
                                            false,
                                            ctx.restriction_,
                                            ctx.time_limit_,
                                            ctx.time_limit_per_supp_);
              int64_t varied_supp = estimated_support - origin_supp;
              if (no_need_to_check) {
                pass_count++;
                //std::cout << "pass gar " << std::endl;
              }
              if (!no_need_to_check || !ctx.support_estimation_) {
                varied_supp = IncHasMatch(pattern,
                                          current_literals,
                                          root_literals[0],
                                          origin_data_graph,
                                          origin_data_graph_nec,
                                          updated_data_graph,
                                          updated_data_graph_nec,
                                          ctx.normal_to_ml_edge_label_,
                                          origin_candidate_set,
                                          origin_candidate_set_removed_nec,
                                          updated_candidate_set,
                                          updated_candidate_set_removed_nec,
                                          id_to_origin_updated_vertex_handle,
                                          updated_vertex_id,
                                          pivot_file_name,
                                          ctx.support_bound_,
                                          false,
                                          ctx.restriction_,
                                          ctx.time_limit_,
                                          ctx.time_limit_per_supp_);
              }
              //std::cout << "varied supp " << varied_supp << std::endl;
              //auto origin_supp = current_literal_handle
              //                    ->template const_attribute<int>(_gar_discover::kSuppKey);
              //std::cout << "origin supp "<< origin_supp << std::endl;
              int64_t updated_supp = origin_supp + varied_supp;
              auto [inc_supp_handle, inc_supp_ret]
                    = current_literal_handle->AddAttribute(kNewSuppKey, origin_supp + varied_supp);
              assert(inc_supp_handle);
              assert(inc_supp_ret);
              if (origin_supp >= ctx.support_bound_
                      && origin_supp + varied_supp < ctx.support_bound_) {
                if (current_literal_handle->CountInEdge() == 0) {
                  continue;
                }
                assert(current_literal_handle->CountInEdge() == 1);
                auto in_edge_it = current_literal_handle->InEdgeBegin();
                auto src_literal_handle = in_edge_it->src_handle();
                child_num_with_valid_supp[src_literal_handle]--;
                if (child_num_with_valid_supp[src_literal_handle] == 0) {
                  // backtrack the literal tree
                  // literal_queue.clear();
                  while (!literal_queue.empty()) {
                    literal_queue.pop();
                  }
                  literal_queue.emplace(root_handle);
                  continue;
                }
              } else if (origin_supp < ctx.support_bound_
                          && origin_supp + varied_supp >= ctx.support_bound_) {
                omp_set_lock(&gar_id_lock);
                auto original_handle = current_literal_handle->FindAttribute(_gar_discover::kOriginalGarID);
                if (original_handle) {
                  auto origin_gar_id = current_literal_handle
                                  ->template const_attribute<int>(_gar_discover::kOriginalGarID);
                  gar_to_extend.push_back(origin_gar_id);
                } else {
                  std::cout << "no original id found " << std::endl;

                }
                omp_unset_lock(&gar_id_lock);
                if (current_literal_handle->CountInEdge() == 0) {
                  continue;
                }
                assert(current_literal_handle->CountInEdge() == 1);
                auto in_edge_it = current_literal_handle->InEdgeBegin();
                auto src_literal_handle = in_edge_it->src_handle();
                child_num_with_valid_supp[src_literal_handle]++;
              }
              if (origin_supp + varied_supp < ctx.support_bound_) {
                // omit changes for all subsequent gars
                continue;
              }
            }

            for (auto out_edge_it = current_literal_handle->OutEdgeBegin();
                    !out_edge_it.IsDone();
                      out_edge_it++) {
              auto dst_literal_handle = out_edge_it->dst_handle();
              literal_queue.emplace(dst_literal_handle);

              auto existed_literals = current_literals;
              DataGraphLiteralType x_literal 
                  = dst_literal_handle->template const_attribute<DataGraphLiteralType>(_gar_discover::kLiteralKey);
              existed_literals.emplace_back(x_literal);
              literal_tree_handle2literals[dst_literal_handle] = existed_literals;
            }
          }
        }
      } else {
        for (auto literal_tree_it = current_generate_tree_node.LiteralTreesBegin();
                  literal_tree_it != current_generate_tree_node.LiteralTreesEnd();
                  literal_tree_it++) {
          auto &literal_tree = *literal_tree_it;
          typename GUNDAM::VertexID<LiteralTreeType>::type kLiteralTreeRootId = 0;

          using LiteralVertexHandleType
                    = typename GUNDAM::VertexHandle<LiteralTreeType>::type;
          
          auto root_handle = literal_tree.FindVertex(kLiteralTreeRootId);
          assert(root_handle->CountInVertex() == 0);

          //std::queue<LiteralVertexHandleType> literal_queue;
          std::queue<std::vector<LiteralVertexHandleType>> literal_queue;
          std::map<LiteralVertexHandleType, std::vector<DataGraphLiteralType>> literal_tree_handle2literals;
          DataGraphLiteralType y_literal 
              = root_handle->template const_attribute<DataGraphLiteralType>(_gar_discover::kLiteralKey);

          std::vector<DataGraphLiteralType> root_literals;
          root_literals.emplace_back(y_literal);

          std::vector<LiteralVertexHandleType> root_literal_vec;
          root_literal_vec.emplace_back(root_handle);
          literal_queue.emplace(root_literal_vec);
          literal_tree_handle2literals[root_handle] = root_literals;

          std::map<LiteralVertexHandleType, unsigned> child_num_with_valid_supp;

          while (!literal_queue.empty()) {
            auto current_literal_handle_vec = literal_queue.front();

            literal_queue.pop();

            std::set<LiteralVertexHandleType> literals_to_continue;

            for (unsigned literal_gar_idx = 0;
                      literal_gar_idx < current_literal_handle_vec.size();
                      literal_gar_idx++) {
              auto current_literal_handle = current_literal_handle_vec[literal_gar_idx];
              auto current_literal_new_supp_handle
                  = current_literal_handle_vec[literal_gar_idx]->FindAttribute(kNewSuppKey);
              if (current_literal_new_supp_handle) {
                auto current_literal_supp = current_literal_handle
                                  ->template const_attribute<int>(_gar_discover::kNewSuppKey);
                if (current_literal_supp < ctx.support_bound_) {
                  continue;
                } else {
                  literals_to_continue.insert(current_literal_handle_vec[literal_gar_idx]);
                }
              } else {
                literals_to_continue.insert(current_literal_handle_vec[literal_gar_idx]);
              }
            }
            std::vector<LiteralVertexHandleType>
                      literals_to_verify_vec(literals_to_continue.begin(),
                                              literals_to_continue.end());
            std::vector<unsigned> valid_child_num(literals_to_verify_vec.size(), 0);

            for (unsigned literal_gar_idx = 0;
                          literal_gar_idx < literals_to_verify_vec.size();
                          literal_gar_idx++) {
              auto& current_literal_handle = literals_to_verify_vec[literal_gar_idx];
              for (auto out_edge_it = current_literal_handle->OutEdgeBegin();
                      !out_edge_it.IsDone();
                        out_edge_it++) {
                auto dst_literal_handle = out_edge_it->dst_handle();
                auto dst_supp = dst_literal_handle
                                  ->template const_attribute<int>(_gar_discover::kSuppKey);
                if (dst_literal_handle->FindAttribute(kNewSuppKey)) {
                  dst_supp = dst_literal_handle
                                  ->template const_attribute<int>(_gar_discover::kNewSuppKey);
                }
                if(dst_supp >= ctx.support_bound_) {
                  valid_child_num[literal_gar_idx]++;
                }
              }
              child_num_with_valid_supp[current_literal_handle] = valid_child_num[literal_gar_idx];
            }
            std::vector<std::string> pivot_file_name_vec;
            std::vector<DataGraphLiteralType> literals_to_check_vec;
            for (auto iter = literals_to_verify_vec.begin();
                      iter != literals_to_verify_vec.end();) {;
              auto current_literal_handle = *iter;

              if (!current_literal_handle) {
                std::cout << "not a valid handle" << std::endl;
              }
              auto current_literal_new_supp_handle = current_literal_handle
                        ->FindAttribute(kNewSuppKey);

              if (child_num_with_valid_supp[current_literal_handle] == 0 && (!current_literal_new_supp_handle)) {

                double probability = current_literal_handle
                                  ->template const_attribute<double>(_gar_discover::kProbability);
                using CandidateSetContainer 
                      = std::map<GraphPatternVertexHandleType,
                              std::vector<DataGraphVertexHandleType>>;
                auto &current_literals = literal_tree_handle2literals[current_literal_handle];

                CandidateSetContainer origin_candidate_set,
                                    updated_candidate_set,
                                      origin_candidate_set_removed_nec,
                                    updated_candidate_set_removed_nec;

                bool init_origin_candidate_set_succ = false;
                bool init_updated_candidate_set_succ = false;

                init_origin_candidate_set_succ
                    = GUNDAM::_dp_iso_using_match
                            ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                      pattern, origin_data_graph, origin_candidate_set);

                
                init_updated_candidate_set_succ
                    = GUNDAM::_dp_iso_using_match
                            ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                      pattern, updated_data_graph, updated_candidate_set);
                

                if (!init_origin_candidate_set_succ) {
                  origin_candidate_set.clear();
                }

                if (!init_updated_candidate_set_succ) {
                  updated_candidate_set.clear();
                }

                bool refine_origin_candidate_set_succ = false;
                bool refine_updated_candidate_set_succ = false;


                refine_origin_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::RefineCandidateSet(
                    pattern, origin_data_graph, origin_candidate_set);


                refine_updated_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::RefineCandidateSet(
                    pattern, updated_data_graph, updated_candidate_set);


                if (!refine_origin_candidate_set_succ) {
                  origin_candidate_set.clear();
                }

                if (!refine_updated_candidate_set_succ) {
                  updated_candidate_set.clear();
                }

                origin_candidate_set_removed_nec = origin_candidate_set;
                GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                    DataGraphType>(
                          origin_candidate_set_removed_nec,
                          origin_data_graph_nec);

                updated_candidate_set_removed_nec = updated_candidate_set;
                GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                    DataGraphType>(
                          updated_candidate_set_removed_nec,
                          updated_data_graph_nec);

                std::string pivot_file_name;
                
                if (ctx.with_index_) {
                  auto pivot_file_handle = current_literal_handle
                            ->FindAttribute(_gar_discover::kMatchFileName);
                  if (!pivot_file_handle) {
                    std::cout << "can not find pivot file name " << std::endl;
                  }
                  pivot_file_name = current_literal_handle
                            ->template const_attribute<std::string>(_gar_discover::kMatchFileName);
                } else {
                  std::cout << "should have indices" << std::endl;
                }
                auto origin_supp = current_literal_handle
                                    ->template const_attribute<int64_t>(_gar_discover::kSuppKey);
                int64_t estimated_support = origin_supp;
                bool no_need_to_check = IncHasMatchByEstimate(pattern,
                                              current_literals,
                                              root_literals[0],
                                              origin_data_graph,
                                              origin_data_graph_nec,
                                              updated_data_graph,
                                              updated_data_graph_nec,
                                              ctx.normal_to_ml_edge_label_,
                                              origin_candidate_set,
                                              origin_candidate_set_removed_nec,
                                              updated_candidate_set,
                                              updated_candidate_set_removed_nec,
                                              id_to_origin_updated_vertex_handle,
                                              updated_vertex_id,
                                              pivot_file_name,
                                              origin_supp,
                                              ctx.support_bound_,
                                              probability,
                                              ctx.threshold_,
                                              estimated_support,
                                              false,
                                              ctx.restriction_,
                                              ctx.time_limit_,
                                              ctx.time_limit_per_supp_);
                if (no_need_to_check && ctx.support_estimation_) {
                  literals_to_verify_vec.erase(iter);
                  auto [inc_supp_handle, inc_supp_ret]
                      = current_literal_handle->AddAttribute(kNewSuppKey, estimated_support);

                  if (origin_supp >= ctx.support_bound_
                          && estimated_support < ctx.support_bound_) {
                    if (current_literal_handle->CountInEdge() == 0) {
                      continue;
                    }
                    assert(current_literal_handle->CountInEdge() == 1);
                    auto in_edge_it = current_literal_handle->InEdgeBegin();
                    auto src_literal_handle = in_edge_it->src_handle();
                    child_num_with_valid_supp[src_literal_handle]--;
                    if (child_num_with_valid_supp[src_literal_handle] == 0) {
                      // backtrack the literal tree
                      // literal_queue.clear();
                      while (!literal_queue.empty()) {
                        literal_queue.pop();
                      }
                      literal_queue.emplace(root_literal_vec);
                      continue;
                    }
                  } else if (origin_supp < ctx.support_bound_
                              && estimated_support >= ctx.support_bound_) {
                    omp_set_lock(&gar_id_lock);
                    auto original_handle = current_literal_handle->FindAttribute(_gar_discover::kOriginalGarID);
                    if (original_handle) {
                      auto origin_gar_id = current_literal_handle
                                      ->template const_attribute<int>(_gar_discover::kOriginalGarID);
                      gar_to_extend.push_back(origin_gar_id);
                    } else {
                      std::cout << "no original id found " << std::endl;
                    }
                    omp_unset_lock(&gar_id_lock);
                    if (current_literal_handle->CountInEdge() == 0) {
                      continue;
                    }

                    assert(current_literal_handle->CountInEdge() == 1);
                    auto in_edge_it = current_literal_handle->InEdgeBegin();
                    auto src_literal_handle = in_edge_it->src_handle();
                    child_num_with_valid_supp[src_literal_handle]++;
                  }
                  if (estimated_support < ctx.support_bound_) {
                    // omit changes for all subsequent gars
                    continue;
                  }
                } else {
                  auto literal_info = (*iter)->template const_attribute<DataGraphLiteralType>(
                                                        _gar_discover::kLiteralKey);
                  literals_to_check_vec.emplace_back(literal_info);
                  pivot_file_name_vec.emplace_back(pivot_file_name);
                  iter++;
                }
              } else {
                iter++;
              }
            }

            if (literals_to_verify_vec.size() == 1) {
                auto& current_literal_handle = literals_to_verify_vec[0];
                auto& current_literals = literal_tree_handle2literals[current_literal_handle];

                using CandidateSetContainer 
                      = std::map<GraphPatternVertexHandleType,
                              std::vector<DataGraphVertexHandleType>>;

                CandidateSetContainer origin_candidate_set,
                                    updated_candidate_set,
                                      origin_candidate_set_removed_nec,
                                    updated_candidate_set_removed_nec;

                bool init_origin_candidate_set_succ = false;
                bool init_updated_candidate_set_succ = false;

                init_origin_candidate_set_succ
                    = GUNDAM::_dp_iso_using_match
                            ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                      pattern, origin_data_graph, origin_candidate_set);

                
                init_updated_candidate_set_succ
                    = GUNDAM::_dp_iso_using_match
                            ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                      pattern, updated_data_graph, updated_candidate_set);
                

                if (!init_origin_candidate_set_succ) {
                  origin_candidate_set.clear();
                }

                if (!init_updated_candidate_set_succ) {
                  updated_candidate_set.clear();
                }

                bool refine_origin_candidate_set_succ = false;
                bool refine_updated_candidate_set_succ = false;


                refine_origin_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::RefineCandidateSet(
                    pattern, origin_data_graph, origin_candidate_set);


                refine_updated_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::RefineCandidateSet(
                    pattern, updated_data_graph, updated_candidate_set);


                if (!refine_origin_candidate_set_succ) {
                  origin_candidate_set.clear();
                }

                if (!refine_updated_candidate_set_succ) {
                  updated_candidate_set.clear();
                }

                origin_candidate_set_removed_nec = origin_candidate_set;
                GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                    DataGraphType>(
                          origin_candidate_set_removed_nec,
                          origin_data_graph_nec);

                updated_candidate_set_removed_nec = updated_candidate_set;
                GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                    DataGraphType>(
                          updated_candidate_set_removed_nec,
                          updated_data_graph_nec);

                std::string pivot_file_name;
                
                if (ctx.with_index_) {
                  auto pivot_file_handle = current_literal_handle
                            ->FindAttribute(_gar_discover::kMatchFileName);
                  if (!pivot_file_handle) {
                    std::cout << "can not find pivot file name " << std::endl;
                  }
                  pivot_file_name = current_literal_handle
                            ->template const_attribute<std::string>(_gar_discover::kMatchFileName);
                }

                auto origin_supp = current_literal_handle
                                    ->template const_attribute<int64_t>(_gar_discover::kSuppKey);

                auto varied_supp = IncHasMatch(pattern,
                                          current_literals,
                                          root_literals[0],
                                          origin_data_graph,
                                          origin_data_graph_nec,
                                          updated_data_graph,
                                          updated_data_graph_nec,
                                          ctx.normal_to_ml_edge_label_,
                                          origin_candidate_set,
                                          origin_candidate_set_removed_nec,
                                          updated_candidate_set,
                                          updated_candidate_set_removed_nec,
                                          id_to_origin_updated_vertex_handle,
                                          updated_vertex_id,
                                          pivot_file_name,
                                          ctx.support_bound_,
                                          false,
                                          ctx.restriction_,
                                          ctx.time_limit_,
                                          ctx.time_limit_per_supp_);
                
              int64_t updated_supp = origin_supp + varied_supp;
              auto [inc_supp_handle, inc_supp_ret]
                    = current_literal_handle->AddAttribute(kNewSuppKey, origin_supp + varied_supp);
              assert(inc_supp_handle);
              assert(inc_supp_ret);
              if (origin_supp >= ctx.support_bound_
                      && origin_supp + varied_supp < ctx.support_bound_) {
                if (current_literal_handle->CountInEdge() == 0) {
                  continue;
                }
                assert(current_literal_handle->CountInEdge() == 1);
                auto in_edge_it = current_literal_handle->InEdgeBegin();
                auto src_literal_handle = in_edge_it->src_handle();
                child_num_with_valid_supp[src_literal_handle]--;
                if (child_num_with_valid_supp[src_literal_handle] == 0) {
                  // backtrack the literal tree
                  // literal_queue.clear();
                  while (!literal_queue.empty()) {
                    literal_queue.pop();
                  }
                  literal_queue.emplace(root_literal_vec);
                  continue;
                }
              } else if (origin_supp < ctx.support_bound_
                          && origin_supp + varied_supp >= ctx.support_bound_) {
                omp_set_lock(&gar_id_lock);
                auto original_handle = current_literal_handle->FindAttribute(_gar_discover::kOriginalGarID);
                if (original_handle) {
                  auto origin_gar_id = current_literal_handle
                                  ->template const_attribute<int>(_gar_discover::kOriginalGarID);
                  gar_to_extend.push_back(origin_gar_id);
                } else {
                  std::cout << "no original id found " << std::endl;
                }
                omp_unset_lock(&gar_id_lock);
                if (current_literal_handle->CountInEdge() == 0) {
                  continue;
                }
                assert(current_literal_handle->CountInEdge() == 1);
                auto in_edge_it = current_literal_handle->InEdgeBegin();
                auto src_literal_handle = in_edge_it->src_handle();
                child_num_with_valid_supp[src_literal_handle]++;
              }
              if (origin_supp + varied_supp < ctx.support_bound_) {
                // omit changes for all subsequent gars
                continue;
              }
            } else if (literals_to_verify_vec.size() > 1) {
              std::vector<DataGraphLiteralType> literals_to_check, common_literals;
              common_literals = root_literals;
              for (auto& literal_handle : literals_to_verify_vec) {
                DataGraphLiteralType x_literal 
                    = literal_handle->template const_attribute<DataGraphLiteralType>(
                                                              _gar_discover::kLiteralKey);
                literals_to_check.emplace_back(x_literal);
              }
              auto in_edge_it = literals_to_verify_vec[0]->InEdgeBegin();
              auto parent_literal_handle = in_edge_it->src_handle();
              common_literals = literal_tree_handle2literals[
                                  (LiteralVertexHandleType)parent_literal_handle];

              using CandidateSetContainer 
                      = std::map<GraphPatternVertexHandleType,
                              std::vector<DataGraphVertexHandleType>>;


                CandidateSetContainer origin_candidate_set,
                                    updated_candidate_set,
                                      origin_candidate_set_removed_nec,
                                    updated_candidate_set_removed_nec;

                bool init_origin_candidate_set_succ = false;
                bool init_updated_candidate_set_succ = false;

                init_origin_candidate_set_succ
                    = GUNDAM::_dp_iso_using_match
                            ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                      pattern, origin_data_graph, origin_candidate_set);

                
                init_updated_candidate_set_succ
                    = GUNDAM::_dp_iso_using_match
                            ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                      pattern, updated_data_graph, updated_candidate_set);
                

                if (!init_origin_candidate_set_succ) {
                  origin_candidate_set.clear();
                }

                if (!init_updated_candidate_set_succ) {
                  updated_candidate_set.clear();
                }

                bool refine_origin_candidate_set_succ = false;
                bool refine_updated_candidate_set_succ = false;


                refine_origin_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::RefineCandidateSet(
                    pattern, origin_data_graph, origin_candidate_set);


                refine_updated_candidate_set_succ
                  = GUNDAM::_dp_iso_using_match
                          ::RefineCandidateSet(
                    pattern, updated_data_graph, updated_candidate_set);


                if (!refine_origin_candidate_set_succ) {
                  origin_candidate_set.clear();
                }

                if (!refine_updated_candidate_set_succ) {
                  updated_candidate_set.clear();
                }

                origin_candidate_set_removed_nec = origin_candidate_set;
                GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                    DataGraphType>(
                          origin_candidate_set_removed_nec,
                          origin_data_graph_nec);

                updated_candidate_set_removed_nec = updated_candidate_set;
                GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                                    DataGraphType>(
                          updated_candidate_set_removed_nec,
                          updated_data_graph_nec);

                std::vector<int64_t> origin_supp_vec;
                for (auto &literal_handle : literals_to_verify_vec) {
                  auto origin_supp = literal_handle
                                    ->template const_attribute<int64_t>(_gar_discover::kSuppKey);
                  origin_supp_vec.push_back(origin_supp);
                }

                auto varied_supp_vec = IncHasMatchForSet(pattern,
                                          common_literals,
                                          root_literals[0],
                                          literals_to_check_vec,
                                          origin_data_graph,
                                          origin_data_graph_nec,
                                          updated_data_graph,
                                          updated_data_graph_nec,
                                          ctx.normal_to_ml_edge_label_,
                                          origin_candidate_set,
                                          origin_candidate_set_removed_nec,
                                          updated_candidate_set,
                                          updated_candidate_set_removed_nec,
                                          id_to_origin_updated_vertex_handle,
                                          updated_vertex_id,
                                          pivot_file_name_vec,
                                          ctx.support_bound_,
                                          false,
                                          ctx.restriction_,
                                          ctx.time_limit_,
                                          ctx.time_limit_per_supp_);

              for (unsigned literal_idx = 0; literal_idx < varied_supp_vec.size(); literal_idx++) {
                int64_t origin_supp = origin_supp_vec[literal_idx];
                int64_t varied_supp = varied_supp_vec[literal_idx];
                auto current_literal_handle = literals_to_verify_vec[literal_idx];
                int64_t updated_supp = origin_supp + varied_supp;
                auto [inc_supp_handle, inc_supp_ret]
                      = current_literal_handle->AddAttribute(kNewSuppKey, origin_supp + varied_supp);
                assert(inc_supp_handle);
                assert(inc_supp_ret);

                if (origin_supp >= ctx.support_bound_
                        && origin_supp + varied_supp < ctx.support_bound_) {
                  if (current_literal_handle->CountInEdge() == 0) {
                    continue;
                  }
                  assert(current_literal_handle->CountInEdge() == 1);
                  auto in_edge_it = current_literal_handle->InEdgeBegin();
                  auto src_literal_handle = in_edge_it->src_handle();
                  child_num_with_valid_supp[src_literal_handle]--;
                  if (child_num_with_valid_supp[src_literal_handle] == 0) {
                    // backtrack the literal tree
                    // literal_queue.clear();
                    while (!literal_queue.empty()) {
                      literal_queue.pop();
                    }
                    literal_queue.emplace(root_literal_vec);
                    continue;
                  }
                } else if (origin_supp < ctx.support_bound_
                            && origin_supp + varied_supp >= ctx.support_bound_) {
                  omp_set_lock(&gar_id_lock);
                  auto original_handle = current_literal_handle->FindAttribute(_gar_discover::kOriginalGarID);
                  if (original_handle) {
                    auto origin_gar_id = current_literal_handle
                                    ->template const_attribute<int>(_gar_discover::kOriginalGarID);
                    gar_to_extend.push_back(origin_gar_id);
                  } else {
                    std::cout << "no original id found " << std::endl;
                  }
                  omp_unset_lock(&gar_id_lock);
                  if (current_literal_handle->CountInEdge() == 0) {
                    continue;
                  }

                  assert(current_literal_handle->CountInEdge() == 1);
                  auto in_edge_it = current_literal_handle->InEdgeBegin();
                  auto src_literal_handle = in_edge_it->src_handle();
                  child_num_with_valid_supp[src_literal_handle]++;
                }
                if (origin_supp + varied_supp < ctx.support_bound_) {
                  // omit changes for all subsequent gars
                  continue;
                }
              }
            }

            for (auto current_literal_handle : literals_to_continue) {
              std::vector<LiteralVertexHandleType> vec_to_push;
              for (auto out_edge_it = current_literal_handle->OutEdgeBegin();
                      !out_edge_it.IsDone();
                        out_edge_it++) {
                auto dst_literal_handle = out_edge_it->dst_handle();
                vec_to_push.emplace_back(dst_literal_handle);
                auto &current_literals = literal_tree_handle2literals[dst_literal_handle];
                auto existed_literals = current_literals;
                DataGraphLiteralType x_literal 
                    = dst_literal_handle->template const_attribute<DataGraphLiteralType>(_gar_discover::kLiteralKey);
                existed_literals.emplace_back(x_literal);
                literal_tree_handle2literals[dst_literal_handle] = existed_literals;
              }
              literal_queue.emplace(vec_to_push);
            }
          }
        }
      }
      if (generate_node_idx % 1000 == 0) {
        std::cout << "pattern generate_node_idx finished" << generate_node_idx << std::endl;
      }
    }


    ExportGeneratedGarsWithSuppAndTree<GraphPatternType, 
                           DataGraphType>(
          0, 
          ctx.generate_tree_level_,
          ctx.restriction_,
          ctx.support_bound_,
          ctx.output_gar_dir_,
          ctx.fid_);

    std::string msg(kInfoProcessPrefix);
    for (unsigned i = 0 ; i < gar_to_extend.size(); i++) {
      msg = msg + " " + std::to_string(gar_to_extend[i]);
    }
    auto& channel_0 = messages.Channels()[0];
    channel_0.SendToFragment(kExpandPatternFragID, msg);
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
    auto& channels = messages.Channels();
    std::vector<int> delivered_gar_ids;
    bool receive_info_process = false,
         receive_deliver_pattern = false,
         receive_message = false;
    messages.ParallelProcess<std::string>(
        // thread_num(),
        1, [&ctx, 
            &delivered_gar_ids,
            &receive_info_process,
            &receive_deliver_pattern,
            &receive_message](int tid, std::string msg) {

          auto res_deliver_pattern =
              std::mismatch(kDeliverPatternPrefix.begin(),
                            kDeliverPatternPrefix. end (), msg.begin());
          auto res_info =     
              std::mismatch(kInfoProcessPrefix.begin(),
                            kInfoProcessPrefix. end (), msg.begin());
          receive_message = true;
          if (res_deliver_pattern.first == kDeliverPatternPrefix.end()) {
            msg = msg.substr(kDeliverPatternPrefix.size());
            // kDeliverPatternPrefix is the prefix of msg.
            receive_deliver_pattern = true;
            // #####################################################################
            // ##  receive the patterns delievered from the kExpandPatternFragID  ##
            // ##  to evaluate the match count                                    ##
            // #####################################################################
              //std::tuple<GraphPatternType, size_t, bool> temp_delivered_pattern;
            if (msg != "") {
              std::stringstream ss(msg);
                int gar_id = 0;
              while (ss >> gar_id) {
                delivered_gar_ids.emplace_back(gar_id);
              }
            }
          } else if (res_info.first == kInfoProcessPrefix.end()) {
            msg = msg.substr(kInfoProcessPrefix.size());
            receive_info_process = true;
            assert(ctx.fid_ == kExpandPatternFragID);
            // ###############################################################
            // ##  receive the number of processors of the ordinary workers  #
            // ###############################################################
            if (msg != "") {
              std::stringstream ss(msg);
              int gar_id = 0;
              while (ss >> gar_id) {
                delivered_gar_ids.emplace_back(gar_id);
              }
            }
          } else {
            // unknown message type
            assert(false);
          }
        });

    auto& origin_data_graph = ctx.data_graph_[0].origin_data_graph();
    auto& updated_data_graph = ctx.data_graph_[0].updated_data_graph();
    auto& origin_data_graph_nec = ctx.data_graph_[0].origin_data_graph_nec();
    auto& updated_data_graph_nec = ctx.data_graph_[0].updated_data_graph_nec();
    auto& id_to_origin_updated_vertex_handle = ctx.data_graph_[0].id_to_vertex_handle();
    auto& updated_vertex_id = ctx.data_graph_[0].updated_vertex_id();
    using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
    using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;
    
    using CandidateSetContainer 
      = std::map<GraphPatternVertexHandleType,
        std::vector<DataGraphVertexHandleType>>;

    if (receive_deliver_pattern) {
      int local_pattern_offset = 0;
      std::map<int, std::vector<DataGraphLiteralType>> pattern_id2lhs_literals;
      std::map<int, std::vector<DataGraphLiteralType>> pattern_id2rhs_literals;
      std::vector<GraphPatternType> pattern_id2pattern;
      std::map<int, int> pattern_id2local_pattern_id;

      int current_worker_id, current_worker_num;
      MPI_Comm_rank(MPI_COMM_WORLD, &current_worker_id);
      MPI_Comm_size(MPI_COMM_WORLD, &current_worker_num);
      omp_lock_t gar_output, expand_pattern_lock;
      omp_init_lock(&gar_output);
      omp_init_lock(&expand_pattern_lock);

      omp_set_nested(false);
      int local_output_count = 0;
      //#pragma omp parallel for schedule(dynamic) 
      for (unsigned i = 0; i < delivered_gar_ids.size(); i++) {
        int gar_id = delivered_gar_ids[i];
        std::vector<DataGraphLiteralType> eligible_literals;
        
        GarType gar_to_extend = ctx.PrepareGARandLiteralsToAdd(gar_id, eligible_literals);

        CandidateSetContainer candidate_set, candidate_set_removed_nec;
        bool init_candidate_set_succ = false;
        init_candidate_set_succ
            = GUNDAM::_dp_iso_using_match
                    ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
              gar_to_extend.pattern(),
              updated_data_graph,
              candidate_set);
        if (!init_candidate_set_succ) {
          continue;
        }

        const bool kRefineCandidateSetSucc 
          = GUNDAM::_dp_iso_using_match
                  ::RefineCandidateSet(
            gar_to_extend.pattern(),
            updated_data_graph, 
            candidate_set);
        if (!kRefineCandidateSetSucc) {
          continue;
        }
        candidate_set_removed_nec = candidate_set;
        GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                            DataGraphType>(
                  candidate_set_removed_nec,
                  updated_data_graph_nec);


        std::vector<DataGraphLiteralType> common_literals;
        for (auto &literal : gar_to_extend.y_literal_set()) {
          common_literals.emplace_back(literal->info());
        }
        for (auto &literal : gar_to_extend.x_literal_set()) {
          common_literals.emplace_back(literal->info());
        }
        std::stack<std::set<int>> literal_stack;

        literal_stack.emplace(std::set<int>());
        while (!literal_stack.empty()) {
          auto existed_literals = literal_stack.top();
          literal_stack.pop();

          std::vector<DataGraphLiteralType> literals_to_check;
          std::vector<int> literals_to_check_id;
          std::vector<DataGraphLiteralType> current_common_literal = common_literals;
          for (unsigned i = 0; i < eligible_literals.size(); i++) {
            if (existed_literals.find(i) != existed_literals.end()) {
              current_common_literal.emplace_back(eligible_literals[i]);
            } else {
              literals_to_check.emplace_back(eligible_literals[i]);
              literals_to_check_id.push_back(i);
            }
          }

          if (current_common_literal.size() >= ctx.j_) {
            continue;
          }

          if (literals_to_check.size() == 0) {
            continue;
          }

          auto [temp_supp_vec, temp_conf_vec]
            = HasMatchForSet(gar_to_extend.pattern(),
                        current_common_literal,
                        *current_common_literal.begin(),
                        literals_to_check,
                        updated_data_graph,
                        ctx.normal_to_ml_edge_label_,
                        candidate_set,
                        candidate_set_removed_nec,
                        ctx.support_bound_, false,
                        ctx.updated_graph_basic_statistics_,
                        ctx.restriction_,
                        ctx.time_limit_,
                        ctx.time_limit_per_supp_);
          

          for (unsigned i = 0; i < literals_to_check.size(); i++) {
            if (temp_supp_vec[i] >= ctx.support_bound_) {
              existed_literals.insert(literals_to_check_id[i]);
              literal_stack.emplace(existed_literals);
              GarType gar_to_output = gar_to_extend;
              for (auto &literal_id : existed_literals) {
                gar_to_output.AddX(eligible_literals[literal_id]);
              }
              omp_set_lock(&gar_output);
              ExportGAR(gar_to_output, ctx.output_gar_dir_, current_worker_id, local_output_count);
              local_output_count++;
              omp_unset_lock(&gar_output);
              existed_literals.erase(literals_to_check_id[i]);
            }
          }
        }

        if (ctx.EligibleForExpandPattern(gar_id)) {
          std::cout << "here " << gar_id << std::endl;
          std::vector<GraphPatternType> expand_pattern = ExpandPattern(gar_to_extend.pattern(),
                      ctx.updated_graph_basic_statistics_, GUNDAM::MaxEdgeId(gar_to_extend.pattern()) + 1,
                      ctx.restriction_, false);
          if (expand_pattern.size() == 0) continue;
          std::vector<DataGraphLiteralType> lhs_literals, rhs_literals;
          ctx.CollectLiterals(gar_id, lhs_literals, rhs_literals);

          omp_set_lock(&expand_pattern_lock);
          for (unsigned pattern_idx = 0; pattern_idx < expand_pattern.size(); pattern_idx++) {
            pattern_id2lhs_literals[local_pattern_offset] = lhs_literals;
            pattern_id2rhs_literals[local_pattern_offset] = rhs_literals;
            pattern_id2pattern.emplace_back(expand_pattern[pattern_idx]);
            pattern_id2local_pattern_id[pattern_idx] = local_pattern_offset;
            local_pattern_offset++;
          }

          omp_unset_lock(&expand_pattern_lock);
        }
      }

      while (pattern_id2pattern.size() != 0) {
        std::vector<GraphPatternType> new_pattern_id2pattern;
        std::map<int, int> new_pattern_id2local_pattern_id;
        for (unsigned pattern_idx = 0; pattern_idx < pattern_id2pattern.size(); pattern_idx++) {

          unsigned local_pattern_idx = pattern_id2local_pattern_id[pattern_idx];
          GarType new_gar = pattern_id2pattern[pattern_idx];
          CandidateSetContainer candidate_set, candidate_set_removed_nec;
          bool init_candidate_set_succ = false;
          init_candidate_set_succ
              = GUNDAM::_dp_iso_using_match
                      ::InitCandidateSet<GUNDAM::MatchSemantics::kIsomorphism>(
                new_gar.pattern(),
                updated_data_graph,
                candidate_set);
          if (!init_candidate_set_succ) {
            continue;
          }

          const bool kRefineCandidateSetSucc 
            = GUNDAM::_dp_iso_using_match
                    ::RefineCandidateSet(
              new_gar.pattern(),
              updated_data_graph, 
              candidate_set);
          if (!kRefineCandidateSetSucc) {
            continue;
          }
          candidate_set_removed_nec = candidate_set;
          GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                              DataGraphType>(
                    candidate_set_removed_nec,
                    updated_data_graph_nec);
          bool eligible_for_expand_pattern = false;
          
          for (auto &rhs_literal : pattern_id2rhs_literals[local_pattern_idx]) {
            new_gar.AddY(rhs_literal);
            std::vector<DataGraphLiteralType> existed_literals;
            existed_literals.emplace_back(rhs_literal);
            auto [temp_supp, temp_conf]
              = HasMatch(new_gar.pattern(),
                          existed_literals,
                        *existed_literals.begin(),
                          updated_data_graph,
                          ctx.normal_to_ml_edge_label_,
                          candidate_set,
                          candidate_set_removed_nec,
                          ctx.support_bound_, false,
                          ctx.updated_graph_basic_statistics_,
                          ctx.restriction_,
                          ctx.time_limit_,
                          ctx.time_limit_per_supp_);
            if (temp_supp < ctx.support_bound_) {
              continue;
            } else {
              eligible_for_expand_pattern = true;
              omp_set_lock(&gar_output);
              ExportGAR(new_gar, ctx.output_gar_dir_, current_worker_id, local_output_count);
              local_output_count++;
              omp_unset_lock(&gar_output);

              if (ctx.j_ > 1) {
                std::stack<std::set<int>> literal_stack;
                literal_stack.emplace(std::set<int>());
                std::vector<DataGraphLiteralType> common_literals;
                common_literals.emplace_back(rhs_literal);

                while (!literal_stack.empty()) {
                  auto existed_literals = literal_stack.top();
                  literal_stack.pop();

                  std::vector<DataGraphLiteralType> literals_to_check;
                  std::vector<int> literals_to_check_id;
                  std::vector<DataGraphLiteralType> current_common_literal = common_literals;
                  for (unsigned i = 0; i < pattern_id2rhs_literals[local_pattern_idx].size(); i++) {
                    if (existed_literals.find(i) != existed_literals.end()) {
                      current_common_literal.emplace_back(pattern_id2rhs_literals[local_pattern_idx][i]);
                    } else {
                      literals_to_check.emplace_back(pattern_id2rhs_literals[local_pattern_idx][i]);
                      literals_to_check_id.push_back(i);
                    }
                  }

                  if (current_common_literal.size() >= ctx.j_) {
                    continue;
                  }

                  if (literals_to_check.size() == 0) {
                    continue;
                  }

                  auto [temp_supp_vec, temp_conf_vec]
                    = HasMatchForSet(new_gar.pattern(),
                                current_common_literal,
                                *current_common_literal.begin(),
                                literals_to_check,
                                updated_data_graph,
                                ctx.normal_to_ml_edge_label_,
                                candidate_set,
                                candidate_set_removed_nec,
                                ctx.support_bound_, false,
                                ctx.updated_graph_basic_statistics_,
                                ctx.restriction_,
                                ctx.time_limit_,
                                ctx.time_limit_per_supp_);
                  

                  for (unsigned i = 0; i < literals_to_check.size(); i++) {
                    if (temp_supp_vec[i] >= ctx.support_bound_) {
                      existed_literals.insert(literals_to_check_id[i]);
                      literal_stack.emplace(existed_literals);
                      GarType gar_to_output = new_gar;
                      for (auto &literal_id : existed_literals) {
                        gar_to_output.AddX(pattern_id2rhs_literals[local_pattern_idx][literal_id]);
                      }
                      omp_set_lock(&gar_output);
                      ExportGAR(new_gar, ctx.output_gar_dir_, current_worker_id, local_output_count);
                      local_output_count++;
                      omp_unset_lock(&gar_output);
                      existed_literals.erase(literals_to_check_id[i]);
                    }
                  }
                }
              }
            }

            if (eligible_for_expand_pattern) {
              std::vector<GraphPatternType> expand_pattern = ExpandPattern(new_gar.pattern(),
                      ctx.updated_graph_basic_statistics_, GUNDAM::MaxEdgeId(new_gar.pattern()) + 1,
                      ctx.restriction_, false);
              if (expand_pattern.size() == 0) continue;

              omp_set_lock(&expand_pattern_lock);
              for (unsigned pattern_idx = 0; pattern_idx < expand_pattern.size(); pattern_idx++) {
                pattern_id2lhs_literals[local_pattern_offset] = pattern_id2lhs_literals[local_pattern_idx];
                pattern_id2rhs_literals[local_pattern_offset] = pattern_id2rhs_literals[local_pattern_idx];
                new_pattern_id2pattern.emplace_back(expand_pattern[pattern_idx]);
                new_pattern_id2local_pattern_id[pattern_idx] = local_pattern_offset;
                local_pattern_offset++;
              }

              omp_unset_lock(&expand_pattern_lock);
            }

          }
        }
        new_pattern_id2pattern.swap(pattern_id2pattern);
        new_pattern_id2local_pattern_id.swap(pattern_id2local_pattern_id);
      }
    } else if (receive_info_process) {
      int current_worker_id, current_worker_num;
      MPI_Comm_rank(MPI_COMM_WORLD, &current_worker_id);
      MPI_Comm_size(MPI_COMM_WORLD, &current_worker_num);
      std::vector<std::vector<int>> gar_ids_to_worker(current_worker_num);

      unsigned w_id = 0;
      for (unsigned i = 0; i < delivered_gar_ids.size(); i++) {
        int gar_id = delivered_gar_ids[i];
        gar_ids_to_worker[w_id].push_back(gar_id);
        w_id++;
        w_id %= current_worker_num;
      }
      std::string msg(kDeliverPatternPrefix);
      std::vector<std::string> msg_to_workers(current_worker_num);
      for (unsigned i = 0; i < gar_ids_to_worker.size(); i++) {
        msg_to_workers[i] = msg;
        for (unsigned j = 0; j < gar_ids_to_worker[i].size(); j++) {
          msg_to_workers[i] = msg_to_workers[i] + " " + std::to_string(gar_ids_to_worker[i][j]);
        }
      }
      auto& channel_0 = messages.Channels()[0];
      for (unsigned dst_fid = 0; dst_fid < current_worker_num; dst_fid++) {
        channel_0.SendToFragment(dst_fid, msg_to_workers[dst_fid]);
      }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> time_span = t2 - ctx.begin_time_;
    std::cout << "total time is " << time_span.count() << std::endl;
    return;
  }
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_INC_GAR_DISCOVER_GAR_DISCOVER_H_
