#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPORT_GAR_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPORT_GAR_H_

#include "gundam/graph_type/small_graph.h"
#include "gundam/graph_type/graph.h"
#include "gundam/graph_type/large_graph.h"
#include "gundam/graph_type/large_graph2.h"

#include "gundam/tool/connected.h"

#include "gar/gar.h"
#include "gar/csv_gar.h"
#include "gar/gpar_new.h"

#include "prob_gar/prob_gar.h"
#include "prob_gar/csv_prob_gar.h"

#include "rule_discover/gar_discover/generate_tree.h"

namespace grape{

namespace _gar_discover {

template<typename GARType>
void ExportGARSet(const std::vector<GARType>& gar_set,
                  const std::string& output_gar_dir,
                  int level,
                  int worker_id,
                  std::string gar_prefix = "gar_level_"){

  const std::string kGarFilePrefix 
            = output_gar_dir 
             + "/" + gar_prefix + std::to_string(level)
                   + "_worker_" + std::to_string(worker_id);

  const std::string v_set_file = kGarFilePrefix + "_v_set.csv";
  const std::string e_set_file = kGarFilePrefix + "_e_set.csv";
  const std::string x_set_file = kGarFilePrefix + "_x_set.csv";
  const std::string y_set_file = kGarFilePrefix + "_y_set.csv";

  gar::WriteGARSet(gar_set, v_set_file, e_set_file,
                            x_set_file, y_set_file);
  return;
}

template<typename GARType>
void ExportGARSetWithTreeAndSupp(const std::vector<GARType>& gar_set,
                  const std::vector<std::pair<int, int>>& tree_supp_set,
                  const std::vector<std::string>& match_file_name_set,
                  const std::vector<double>& probability_set,
                  const std::string& output_gar_dir,
                  int level,
                  int worker_id,
                  std::string gar_prefix = "gar_level_"){

  const std::string kGarFilePrefix 
            = output_gar_dir 
             + "/" + gar_prefix + std::to_string(level)
                   + "_worker_" + std::to_string(worker_id);

  const std::string v_set_file = kGarFilePrefix + "_v_set.csv";
  const std::string e_set_file = kGarFilePrefix + "_e_set.csv";
  const std::string x_set_file = kGarFilePrefix + "_x_set.csv";
  const std::string y_set_file = kGarFilePrefix + "_y_set.csv";

  const std::string tree_supp_set_file = kGarFilePrefix + "_t_set.csv";

  std::ofstream t_output(tree_supp_set_file);
  if (match_file_name_set.empty()) {
    for (unsigned i = 0; i < tree_supp_set.size(); i++) {
      t_output << worker_id << " " << level << " " << i << " "
              << tree_supp_set[i].first << " " << tree_supp_set[i].second
              << "\n";
    }
  } else {
    for (unsigned i = 0; i < tree_supp_set.size(); i++) {
      t_output << worker_id << " " << level << " " << i << " "
              << tree_supp_set[i].first << " " << tree_supp_set[i].second
              << " " << match_file_name_set[i]
              << " " << probability_set[i] << "\n";
    }
  }

  t_output.close();

  gar::WriteGARSet(gar_set, v_set_file, e_set_file,
                            x_set_file, y_set_file);
  return;
}



template<typename ProbGARType>
void ExportProbGARSet(const std::vector<ProbGARType>& prob_gar_set,
                      const std::string& output_gar_dir,
                      int level,
                      int worker_id,
                      std::string gar_prefix = "gar_level_"){

  const std::string kGarFilePrefix 
            = output_gar_dir 
             + "/" + gar_prefix + std::to_string(level)
                   + "_worker_" + std::to_string(worker_id);

  const std::string v_set_file = kGarFilePrefix + "_v_set.csv";
  const std::string e_set_file = kGarFilePrefix + "_e_set.csv";
  const std::string x_set_file = kGarFilePrefix + "_x_set.csv";
  const std::string y_set_file = kGarFilePrefix + "_y_set.csv";
  const std::string p_set_file = kGarFilePrefix + "_p_set.csv";

  prob_gar::WriteProbGARSet(prob_gar_set, v_set_file, e_set_file,
                                          x_set_file, y_set_file,
                                                      p_set_file);
  return;
}

template<typename GraphPatternType,
         typename    DataGraphType,
         typename    SubstructureLevelType>
void ExportGeneratedGars(
      const SubstructureLevelType& current_pattern_edge_size,
      const GenerateTreeLevel<GraphPatternType,
                                 DataGraphType>& current_generate_level,
            const Restriction<GraphPatternType,
                                 DataGraphType>& restriction,
      const std::string& kOutputGarDir,
      const int worker_id){

  using DataGraphTypeLiteralType = gar::LiteralInfo<GraphPatternType, 
                                                       DataGraphType>;

  using GARType = gar::GraphAssociationRule<GraphPatternType, 
                                               DataGraphType>;

  using ProbGarType = prob_gar::ProbGraphAssociationRule<
                                    GraphPatternType, 
                                       DataGraphType>;

  util::Info("export generate level: "
            + std::to_string(current_pattern_edge_size));

  // for not specifed confidence bound
  std::vector<GARType> gar_set;

  // for specifed confidence bound
  std::vector<ProbGarType> prob_gar_set;

  // #pragma omp parallel for schedule(dynamic)
  for (size_t current_generate_node_idx = 0; 
              current_generate_node_idx < static_cast<int>(current_generate_level.size()); 
              current_generate_node_idx++) {

    const auto& generate_node = current_generate_level.node(current_generate_node_idx);

    const auto& generated_pattern = generate_node.const_pattern();

    // assert(GUNDAM::Connected(generated_pattern));

    for (auto literal_tree_it  = generate_node.LiteralTreesBegin();
              literal_tree_it != generate_node.LiteralTreesEnd();
              literal_tree_it++) {

      auto& literal_tree = *literal_tree_it;

      constexpr typename GUNDAM::VertexID<LiteralTreeType>::type kLiteralTreeRootId = 0;

      using LiteralVertexHandleType = typename GUNDAM::VertexHandle<const LiteralTreeType>::type;
      
      auto root_handle = literal_tree.FindVertex(kLiteralTreeRootId);
      assert(root_handle->CountInVertex() == 0);
      std::stack<std::pair<LiteralVertexHandleType, bool>> literal_stack;

      DataGraphTypeLiteralType y_literal 
          = root_handle->template const_attribute<DataGraphTypeLiteralType>(_gar_discover::kLiteralKey);

      GARType gar(generated_pattern);

      gar.AddY(y_literal);

      gar_set.emplace_back(gar);

      if (restriction.specified_confidence_bound()) {
        float confidence
            = root_handle->template const_attribute<float>(_gar_discover::kConfKey);
        if (confidence >= restriction.confidence_bound()) {
          prob_gar_set.emplace_back(gar, confidence);
        }
      }

      #ifndef NDEBUG
      // used to find whether the same vertex has been visited twice
      std::set<LiteralVertexHandleType> visited;
      visited.emplace(root_handle);
      #endif // NDEBUG
      
      for (auto out_vertex_cit = root_handle->OutVertexBegin();
               !out_vertex_cit.IsDone();
                out_vertex_cit++) {
        literal_stack.emplace(out_vertex_cit, true);
        #ifndef NDEBUG
        visited.emplace(out_vertex_cit);
        #endif // NDEBUG
      }

      std::vector<DataGraphTypeLiteralType> existed_literals;

      while (!literal_stack.empty()) {
        auto [ current_literal_handle,
               current_literal_considered_first_time ] 
                     = literal_stack.top();
        literal_stack.pop();
        if (!current_literal_considered_first_time) {
          // recycle
          // the label of the vertex in the literal tree represents
          // the level of it 
          assert(existed_literals.size() 
              == current_literal_handle->label() - 1); // literal held in root does not contained in it
          existed_literals.pop_back();
          continue;
        }

        literal_stack.emplace(current_literal_handle, false);
        // expanding
        existed_literals.emplace_back(current_literal_handle
                ->template const_attribute<
                                DataGraphTypeLiteralType>(_gar_discover::kLiteralKey));

        // the label of the vertex in the literal tree represents
        // the level of it
        assert(existed_literals.size() 
            == current_literal_handle->label() - 1); // literal held in root does not contained in it

        GARType temp_gar(gar);
        // generate a new gar
        for (const auto& literal : existed_literals){
          temp_gar.AddX(literal);
        }
        gar_set.emplace_back(temp_gar);

        if (restriction.specified_confidence_bound()) {
          float confidence
              = current_literal_handle->template const_attribute<float>(_gar_discover::kConfKey);
          if (confidence >= restriction.confidence_bound()) {
            prob_gar_set.emplace_back(temp_gar, confidence);
          }
        }

        for (auto out_literal_vertex_cit = current_literal_handle->OutVertexBegin();
                 !out_literal_vertex_cit.IsDone();
                  out_literal_vertex_cit++) {
          // should not have visited this vertex before
          assert(visited.find(out_literal_vertex_cit) == visited.cend());
          // add to stack
          literal_stack.emplace(out_literal_vertex_cit, true);
        }
      }
    }
  }
  util::Info("gar_set.size(): " + std::to_string(gar_set.size()));
  if (!restriction.specified_confidence_bound()) {
    if (gar_set.empty()) {
      return;
    }
    ExportGARSet(gar_set, kOutputGarDir, current_pattern_edge_size, worker_id);
    return;
  }
  if (prob_gar_set.empty()) {
    return;
  }
  ExportProbGARSet(prob_gar_set, kOutputGarDir, current_pattern_edge_size, worker_id);
  return;
}

template<typename GraphPatternType,
         typename    DataGraphType>
void ExportGAR(const gar::GraphAssociationRule<GraphPatternType, 
                                               DataGraphType>& gar,
              const std::string& kOutputGarDir,
              const int worker_id,
              const int gar_idx) {
  using DataGraphTypeLiteralType = gar::LiteralInfo<GraphPatternType, 
                                                      DataGraphType>;
  std::string output_v_file_address = kOutputGarDir + "/" + "inc_gar_worker_"
              + std::to_string(worker_id) +"_"
              + std::to_string(gar_idx) + "_v.csv";
  std::string output_e_file_address = kOutputGarDir + "/" + "inc_gar_worker_"
              + std::to_string(worker_id) +"_"
              + std::to_string(gar_idx) + "_e.csv";
  std::string output_x_file_address = kOutputGarDir + "/" + "inc_gar_worker_"
              + std::to_string(worker_id) +"_"
              + std::to_string(gar_idx) + "_x.csv";
  std::string output_y_file_address = kOutputGarDir + "/" + "inc_gar_worker_"
              + std::to_string(worker_id) +"_"
              + std::to_string(gar_idx) + "_y.csv";

  gar::WriteGAR(gar, output_v_file_address, output_e_file_address,
                     output_x_file_address, output_y_file_address);
  return;
}


template<typename GraphPatternType,
         typename    DataGraphType,
         typename    SubstructureLevelType>
void ExportGeneratedGarsWithSuppAndTree(
      const SubstructureLevelType& current_pattern_edge_size,
      const GenerateTreeLevel<GraphPatternType,
                                 DataGraphType>& current_generate_level,
            const Restriction<GraphPatternType,
                                 DataGraphType>& restriction,
            const uint64_t supp_bound,
      const std::string& kOutputGarDir,
      const int worker_id){

  using DataGraphTypeLiteralType = gar::LiteralInfo<GraphPatternType, 
                                                       DataGraphType>;

  using GARType = gar::GraphAssociationRule<GraphPatternType, 
                                               DataGraphType>;

  using ProbGarType = prob_gar::ProbGraphAssociationRule<
                                    GraphPatternType, 
                                       DataGraphType>;

  util::Info("export generate level with supp and tree: "
            + std::to_string(current_pattern_edge_size));

  // for not specifed confidence bound
  std::vector<GARType> gar_set;

  std::vector<std::pair<int, int>> tree_supp_set;
  std::vector<std::string> match_file_name_set;
  std::vector<double> probability_set;

  // for specifed confidence bound
  std::vector<ProbGarType> prob_gar_set;

  // #pragma omp parallel for schedule(dynamic)
  for (size_t current_generate_node_idx = 0; 
              current_generate_node_idx < static_cast<int>(current_generate_level.size()); 
              current_generate_node_idx++) {

    const auto& generate_node = current_generate_level.node(current_generate_node_idx);

    const auto& generated_pattern = generate_node.const_pattern();

    // assert(GUNDAM::Connected(generated_pattern));

    for (auto literal_tree_it  = generate_node.LiteralTreesBegin();
              literal_tree_it != generate_node.LiteralTreesEnd();
              literal_tree_it++) {

      auto& literal_tree = *literal_tree_it;

      constexpr typename GUNDAM::VertexID<LiteralTreeType>::type kLiteralTreeRootId = 0;

      using LiteralVertexHandleType = typename GUNDAM::VertexHandle<const LiteralTreeType>::type;
      
      auto root_handle = literal_tree.FindVertex(kLiteralTreeRootId);
      assert(root_handle->CountInVertex() == 0);
      std::stack<std::pair<LiteralVertexHandleType, bool>> literal_stack;
      std::stack<int> tree_offset_stack; // pair<parent_offset, cur_offset>

      DataGraphTypeLiteralType y_literal 
          = root_handle->template const_attribute<DataGraphTypeLiteralType>(_gar_discover::kLiteralKey);

      GARType gar(generated_pattern);

      gar.AddY(y_literal);

      gar_set.emplace_back(gar);
      int supp
            = root_handle->template const_attribute<int>(_gar_discover::kSuppKey);
      tree_supp_set.emplace_back(-1, supp);

      auto root_name_handle = root_handle->FindAttribute(_gar_discover::kMatchFileName);
      if (root_name_handle) {
        match_file_name_set.emplace_back(root_handle->
                        template const_attribute<std::string>(_gar_discover::kMatchFileName));
        //std::cout << "here with string " << match_file_name_set.back() << std::endl;
      }

      auto root_probability_handle = root_handle->FindAttribute(_gar_discover::kProbability);
      if (root_probability_handle) {
        probability_set.emplace_back(root_handle->
                        template const_attribute<double>(_gar_discover::kProbability));
        //std::cout << "here with string " << match_file_name_set.back() << std::endl;
      }

      if (restriction.specified_confidence_bound()) {
        float confidence
            = root_handle->template const_attribute<float>(_gar_discover::kConfKey);
        if (confidence >= restriction.confidence_bound()) {
          prob_gar_set.emplace_back(gar, confidence);
        }
      }

      #ifndef NDEBUG
      // used to find whether the same vertex has been visited twice
      std::set<LiteralVertexHandleType> visited;
      visited.emplace(root_handle);
      #endif // NDEBUG
      
      for (auto out_vertex_cit = root_handle->OutVertexBegin();
               !out_vertex_cit.IsDone();
                out_vertex_cit++) {
        literal_stack.emplace(out_vertex_cit, true);
        tree_offset_stack.push(gar_set.size() - 1);
        #ifndef NDEBUG
        visited.emplace(out_vertex_cit);
        #endif // NDEBUG
      }

      std::vector<DataGraphTypeLiteralType> existed_literals;

      while (!literal_stack.empty()) {
        auto [ current_literal_handle,
               current_literal_considered_first_time ] 
                     = literal_stack.top();
        int current_parent_offset = tree_offset_stack.top();

        literal_stack.pop();
        tree_offset_stack.pop();
        if (!current_literal_considered_first_time) {
          // recycle
          // the label of the vertex in the literal tree represents
          // the level of it 
          assert(existed_literals.size() 
              == current_literal_handle->label() - 1); // literal held in root does not contained in it
          existed_literals.pop_back();
          continue;
        }

        literal_stack.emplace(current_literal_handle, false);
        tree_offset_stack.push(current_parent_offset);
        // expanding
        existed_literals.emplace_back(current_literal_handle
                ->template const_attribute<
                                DataGraphTypeLiteralType>(_gar_discover::kLiteralKey));

        // the label of the vertex in the literal tree represents
        // the level of it
        assert(existed_literals.size() 
            == current_literal_handle->label() - 1); // literal held in root does not contained in it

        GARType temp_gar(gar);
        // generate a new gar
        for (const auto& literal : existed_literals){
          temp_gar.AddX(literal);
        }
        gar_set.emplace_back(temp_gar);

        int curr_supp = current_literal_handle->
                template const_attribute<int>(_gar_discover::kSuppKey);
        auto new_supp_attr_handle = current_literal_handle->FindAttribute(_gar_discover::kNewSuppKey);
        if (new_supp_attr_handle) {
          curr_supp = current_literal_handle
                          ->template const_attribute<int>(_gar_discover::kNewSuppKey);
        }

        tree_supp_set.emplace_back(current_parent_offset, curr_supp);

        auto match_file_name_handle = current_literal_handle->FindAttribute(_gar_discover::kMatchFileName);
        if (match_file_name_handle) {
          match_file_name_set.emplace_back(current_literal_handle
                          ->template const_attribute<std::string>(_gar_discover::kMatchFileName));
          //std::cout << "here with string " << match_file_name_set.back() << std::endl;
        }

        auto probability_handle = current_literal_handle->FindAttribute(_gar_discover::kProbability);
        if (probability_handle) {
          probability_set.emplace_back(current_literal_handle
                          ->template const_attribute<double>(_gar_discover::kProbability));
          //std::cout << "here with string " << match_file_name_set.back() << std::endl;
        }

        if (restriction.specified_confidence_bound()) {
          float confidence
              = current_literal_handle->template const_attribute<float>(_gar_discover::kConfKey);
          if (confidence >= restriction.confidence_bound()) {
            prob_gar_set.emplace_back(temp_gar, confidence);
          }
        }

        if (curr_supp < supp_bound) {
          continue;
        }

        for (auto out_literal_vertex_cit = current_literal_handle->OutVertexBegin();
                 !out_literal_vertex_cit.IsDone();
                  out_literal_vertex_cit++) {
          // should not have visited this vertex before
          assert(visited.find(out_literal_vertex_cit) == visited.cend());
          // add to stack
          literal_stack.emplace(out_literal_vertex_cit, true);
          tree_offset_stack.emplace(gar_set.size() - 1);
        }
      }
    }
  }
  if (current_pattern_edge_size == 0) {
    // for inc gar discover
    match_file_name_set.clear();
  }
  if ((tree_supp_set.size() != match_file_name_set.size())
        && (!match_file_name_set.empty())) {
    std::cout << "file name set size should be same size or empty!" << std::endl;
    return;
  }
  util::Info("gar_set.size(): " + std::to_string(gar_set.size()));
  if (!restriction.specified_confidence_bound()) {
    if (gar_set.empty()) {
      return;
    }
    ExportGARSetWithTreeAndSupp(gar_set, tree_supp_set, match_file_name_set, probability_set,
          kOutputGarDir, current_pattern_edge_size, worker_id);
    return;
  }
//  if (prob_gar_set.empty()) {
//    return;
//  }
  ExportGARSetWithTreeAndSupp(gar_set, tree_supp_set, match_file_name_set, probability_set,
                               kOutputGarDir, current_pattern_edge_size, worker_id);
  return;
}


template<typename      GraphPatternType,
         typename         DataGraphType,
         typename SubstructureLevelType>
std::string ExportGARMatch(const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal_info,
                          std::vector<std::pair<
                                  std::vector<typename GUNDAM::VertexID<DataGraphType>::type>,
                                  std::vector<typename GUNDAM::EdgeID<DataGraphType>::type>>>& match_vec,
                          SubstructureLevelType current_level,
                          uint64_t pattern_id, uint64_t gar_offset, int worker_id,
                          std::string &output_gar_dir) {
  std::string output_file_name = output_gar_dir + "/gar_match_worker_" + std::to_string(worker_id)
                               + "_level_"   + std::to_string(current_level)
                               + "_pattern_" + std::to_string(pattern_id)
                               + "_gar_offset_" + std::to_string(gar_offset);
  //std::cout << "output file name " << output_file_name << std::endl;
  std::ofstream output_file(output_file_name);
  if (!output_file.is_open()) {
    std::cout << "can not open file " << output_file_name << std::endl;
    return output_file_name;
  }
  auto x_id = rhs_literal_info.x_id();
  auto y_id = rhs_literal_info.y_id();
  output_file << x_id << " " << y_id << "\n";
  for (unsigned match_idx = 0; match_idx < match_vec.size(); match_idx++) {
    auto &pivot_vec = match_vec[match_idx].first;
    auto &edge_vec = match_vec[match_idx].second;
    for (auto &pivot_id : pivot_vec) {
      output_file << " " << pivot_id ;
    }
    for (auto &edge_id : edge_vec) {
      output_file << " " << edge_id;
    }
    output_file << "\n";
  }
  output_file.close();
  return output_file_name;
}

} // namespace _gar_discover

} // namespace grape

#endif // _EXPORT_GAR_H
