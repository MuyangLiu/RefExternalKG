#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPAND_LITERAL_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPAND_LITERAL_H_

#include "rule_discover/gar_discover/expand_tree.h"
#include "rule_discover/gar_discover/restriction.h"

#include "gundam/graph_type/graph.h"
#include "gundam/graph_type/large_graph.h"
#include "gundam/graph_type/large_graph2.h"
#include "gundam/graph_type/small_graph.h"

#include "gundam/tool/has_edge.h"
#include "gundam/tool/distance.h"
#include "gundam/tool/connected_component.h"
#include "gundam/tool/topological/star/is_star.h"

#include "gundam/type_getter/vertex_label.h"
#include "gundam/type_getter/edge_label.h"

#include "gundam/graph_statistics/graph_basic_statistics.h"

#include "gar/literal_info.h"

namespace grape {

namespace _gar_discover {

// #################
// #  optimize me  #
// #################
// expand rhs literal
template <typename GraphPatternType,
          typename    DataGraphType>
void _ExpandLiteral(
    const ExpandTreeNode<GraphPatternType,
                            DataGraphType>& expand_tree_node,
    const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& vertex_handle_set, // each literal needs to contain
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
    const std::set<typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_edge_label_set, 
    const std::set<
        std::tuple<typename GUNDAM::VertexLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::VertexLabel<GraphPatternType>::type>>& ml_edge_type_set,
    const std::map<typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_to_normal_edge_label,
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
                                  const bool rhs_literal_only,
          std::vector<
           gar::LiteralInfo<GraphPatternType,
                               DataGraphType>>& expanded_literal_info) {
                                
  using PatternVertexHandleType      = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using PatternVertexIDType          = typename GUNDAM::VertexID    <GraphPatternType>::type;
  using PatternVertexConstHandleType = typename GUNDAM::VertexHandle<const GraphPatternType>::type;

  assert(vertex_handle_set.size() <= 2);

  // all vertexes contained in literals need to be contained
  // in vertex_handle_set
  std::vector<PatternVertexHandleType> 
            sorted_vertex_handle_set(vertex_handle_set);
  std::sort(sorted_vertex_handle_set.begin(),
            sorted_vertex_handle_set.end());
  assert(std::is_sorted(sorted_vertex_handle_set.begin(),
                        sorted_vertex_handle_set.end()));

  // at least one vertex needs to be contained in sorted_vertex_id_set
  std::vector<PatternVertexIDType> sort_new_vertex_id_set; // each literal must contain at least one
  for (auto new_vertex_id_cit  = expand_tree_node.NewVertexesIdCBegin();
            new_vertex_id_cit != expand_tree_node.NewVertexesIdCEnd();
            new_vertex_id_cit++) {
    sort_new_vertex_id_set.emplace_back(*new_vertex_id_cit);
  }
  std::sort(sort_new_vertex_id_set.begin(),
            sort_new_vertex_id_set.end());
  assert(std::is_sorted(sort_new_vertex_id_set.begin(),
                        sort_new_vertex_id_set.end()));

  /* ############### *
   * ##  for GCR  ## *
   * ############### */
  std::vector<
  std::vector<PatternVertexIDType>> components_vertex_id_set;
  std::vector<PatternVertexIDType>  star_end_vertex_id_set;

  if (restriction.gcr()) {
    components_vertex_id_set.resize(2);
    std::vector<GraphPatternType> connected_components 
                        = GUNDAM::ConnectedComponent(expand_tree_node.const_pattern());
    assert(connected_components.size() == 2);
    for (size_t cc_idx = 0; cc_idx < 2 /* at most have two star */; cc_idx++) {
      auto& vertex_id_set = components_vertex_id_set[cc_idx];
      const auto& cc = connected_components[cc_idx];
      assert(GUNDAM::IsStar<true>(cc));
      vertex_id_set.reserve(cc.CountVertex());
      for (auto vertex_it = cc.VertexBegin();
               !vertex_it.IsDone();
                vertex_it++) {
        vertex_id_set.emplace_back(vertex_it->id());
      }
      assert(vertex_id_set.size() == cc.CountVertex());
      std::sort(vertex_id_set.begin(),
                vertex_id_set.end());
      
      // add to star end points 
      const auto [end_vertex_handle_set,
              central_vertex_handle] = GUNDAM::StarEndPoints<true>(cc);
      // should be a star
      assert(end_vertex_handle_set.size() >= 2);
      assert(end_vertex_handle_set.size() == 2 || central_vertex_handle);
      for (const auto& vertex_handle 
                 : end_vertex_handle_set) {
        star_end_vertex_id_set.emplace_back(vertex_handle->id());
      }
    }
    std::sort(star_end_vertex_id_set.begin(),
              star_end_vertex_id_set.end());
  }

  if (restriction.consider_variable_literal()
   &&(sorted_vertex_handle_set.empty() // has not specified vertex set 
   || sorted_vertex_handle_set.size() == 2)) { // has specified two vertexes
    // ###############################
    // ##  variable literal         ##
    // ##  x.A = y.B                ##
    // ##  enumerate x, y, A and B  ##
    // ###############################
    for (auto x_vertex_cit = expand_tree_node.const_pattern()
                                             .VertexBegin(); 
             !x_vertex_cit.IsDone(); 
              x_vertex_cit++) {
      PatternVertexConstHandleType x_vertex_handle = x_vertex_cit;
      // enumerate x vertex
      if (!sorted_vertex_handle_set.empty()
       && !std::binary_search(sorted_vertex_handle_set.begin(),
                              sorted_vertex_handle_set.end(),
                                   x_vertex_handle)) {
        // specified vertex set but does not contain x_vertex_cit
        continue;
      }
      if (restriction.gcr()) {
        if (!std::binary_search(star_end_vertex_id_set.begin(),
                                star_end_vertex_id_set.end(),
                                       x_vertex_handle->id())) {
          // is not an end point
          continue;
        }
      }
      for (auto y_vertex_cit = x_vertex_cit; 
               !y_vertex_cit.IsDone(); 
                y_vertex_cit++) {
        if (x_vertex_cit == y_vertex_cit) {
          // enumerate y vertex
          // does not consider two attributes of the same vertex
          continue;
        }
        PatternVertexConstHandleType y_vertex_handle = y_vertex_cit;
        if (restriction.gcr()) {
          if (!std::binary_search(star_end_vertex_id_set.begin(),
                                  star_end_vertex_id_set.end(),
                                         y_vertex_handle->id())) {
            // is not an end point
            continue;
          }
          // to check whether x and y are in the
          // same connected_component
          bool contained_in_same_cc = false;
          for (size_t cc_idx = 0; cc_idx < 2 /* at most have two star */; 
                      cc_idx++) {
            auto& vertex_id_set = components_vertex_id_set[cc_idx];
            if (std::binary_search(vertex_id_set.begin(),
                                   vertex_id_set.end(),
                                 x_vertex_handle->id())) {
              // to find whether the other literal is also contained
              // in this component
              if (std::binary_search(vertex_id_set.begin(),
                                     vertex_id_set.end(),
                                   y_vertex_handle->id())) {
                assert(!std::binary_search(components_vertex_id_set[1 - cc_idx].begin(),
                                           components_vertex_id_set[1 - cc_idx].end(),
                                                    y_vertex_handle->id()));
                contained_in_same_cc = true;
                break;
              }
              assert(std::binary_search(components_vertex_id_set[1 - cc_idx].begin(),
                                        components_vertex_id_set[1 - cc_idx].end(),
                                                 y_vertex_handle->id()));
              assert(!contained_in_same_cc);
              break;
            }
          }
          if (contained_in_same_cc) {
            continue;
          }
        }
        if (!sorted_vertex_handle_set.empty()
         && !std::binary_search(sorted_vertex_handle_set.begin(),
                                sorted_vertex_handle_set.end(),
                                     y_vertex_handle)) {
          // specified vertex set but does not contain y_vertex_cit
          continue;
        }
        if (!std::binary_search(sort_new_vertex_id_set.begin(),
                                sort_new_vertex_id_set.end(),
                                       x_vertex_cit->id())
         && !std::binary_search(sort_new_vertex_id_set.begin(),
                                sort_new_vertex_id_set.end(),
                                       y_vertex_cit->id())) {
          // both x and y are not contained in the new vertexes
          // move to the next vertex
          continue;
        }
        if (restriction.variable_literal_only_between_connected_vertexes()) {
          if (!GUNDAM::HasEdge<const GraphPatternType, true>(
                                    x_vertex_cit,
                                    y_vertex_cit)) {
            continue;
          }
        }
        if (restriction.variable_literal_only_between_2_hop_connected_vertexes()) {
          auto [distance, connected] = GUNDAM::Distance<true>(expand_tree_node.pattern(),
                                                              x_vertex_cit->id(),
                                                              y_vertex_cit->id());
          assert(connected);
          if (!connected || distance > 2) {
            continue;
          }
        }
        // two attributes at different vertexes
        // enumerate A, e.g. the attribute key of x
         auto  x_attr_ret  = graph_basic_statistics.legal_attr_set().find(x_vertex_cit->label());
        assert(x_attr_ret != graph_basic_statistics.legal_attr_set().end());
        for (const auto& [x_attr_key, x_attr_value_set] : x_attr_ret->second) {
          // enumerate B, e.g. the attribute key of y
           auto  y_attr_ret  = graph_basic_statistics.legal_attr_set().find(y_vertex_cit->label());
          assert(y_attr_ret != graph_basic_statistics.legal_attr_set().end());
          for (const auto& [y_attr_key, y_attr_value_set] : y_attr_ret->second) {
            if (x_attr_key != y_attr_key){
              // only add variable literal with same key
              continue;
            }
            expanded_literal_info.emplace_back(x_vertex_cit->id(), x_attr_key, 
                                               y_vertex_cit->id(), y_attr_key);
          }
        }
      }
    }
  }

  if (restriction.consider_edge_literal()
   &&(sorted_vertex_handle_set.empty() // has not specified vertex set 
   || sorted_vertex_handle_set.size() == 2)) { // has specified two vertexes
    // ###########################
    // ##  edge literal         ##
    // ##  x - edge_label -> y  ##
    // ###########################
    for (auto x_vertex_cit = expand_tree_node.const_pattern()
                                             .VertexBegin(); 
             !x_vertex_cit.IsDone(); 
              x_vertex_cit++) {
      PatternVertexConstHandleType x_vertex_handle = x_vertex_cit;
      if (!sorted_vertex_handle_set.empty()
       && !std::binary_search(sorted_vertex_handle_set.begin(),
                              sorted_vertex_handle_set.end(),
                                   x_vertex_handle)) {
        // specified vertex set but does not contain x_vertex_cit
        continue;
      }
      if (restriction.gcr()) {
        if (!std::binary_search(star_end_vertex_id_set.begin(),
                                star_end_vertex_id_set.end(),
                                       x_vertex_handle->id())) {
          // is not an end point
          continue;
        }
      }
      for (auto y_vertex_cit = expand_tree_node.const_pattern()
                                               .VertexBegin(); 
               !y_vertex_cit.IsDone(); 
                y_vertex_cit++) {
        // not add self loop as edge literal
        if (x_vertex_cit == y_vertex_cit) 
          continue;
          
        PatternVertexConstHandleType y_vertex_handle = y_vertex_cit;

        if (restriction.gcr()) {
          if (!std::binary_search(star_end_vertex_id_set.begin(),
                                  star_end_vertex_id_set.end(),
                                         y_vertex_handle->id())) {
            // is not an end point
            continue;
          }
          // to check whether x and y are in the
          // same connected_component
          bool contained_in_same_cc = false;
          for (size_t cc_idx = 0; cc_idx < 2 /* at most have two star */; 
                      cc_idx++) {
            auto& vertex_id_set = components_vertex_id_set[cc_idx];
            if (std::binary_search(vertex_id_set.begin(),
                                   vertex_id_set.end(),
                                 x_vertex_handle->id())) {
              // to find whether the other literal is also contained
              // in this component
              if (std::binary_search(vertex_id_set.begin(),
                                     vertex_id_set.end(),
                                   y_vertex_handle->id())) {
                assert(!std::binary_search(components_vertex_id_set[1 - cc_idx].begin(),
                                           components_vertex_id_set[1 - cc_idx].end(),
                                                    y_vertex_handle->id()));
                contained_in_same_cc = true;
                break;
              }
              assert(std::binary_search(components_vertex_id_set[1 - cc_idx].begin(),
                                        components_vertex_id_set[1 - cc_idx].end(),
                                                 y_vertex_handle->id()));
              assert(!contained_in_same_cc);
              break;
            }
          }
          if (contained_in_same_cc) {
            continue;
          }
        }

        if (!sorted_vertex_handle_set.empty()
         && !std::binary_search(sorted_vertex_handle_set.begin(),
                                sorted_vertex_handle_set.end(),
                                     y_vertex_handle)) {
          // specified vertex set but does not contain y_vertex_cit
          continue;
        }

        if (!std::binary_search(sort_new_vertex_id_set.begin(),
                                sort_new_vertex_id_set.end(),
                                       x_vertex_cit->id())
         && !std::binary_search(sort_new_vertex_id_set.begin(),
                                sort_new_vertex_id_set.end(),
                                       y_vertex_cit->id())) {
          // both x and y are not contained in the new vertexes
          // move to the next vertex
          continue;
        }
        
        if (restriction.edge_literal_only_between_2_hop_connected_vertexes()) {
          auto [distance, connected] = GUNDAM::Distance<true>(expand_tree_node.pattern(),
                                                              x_vertex_cit->id(),
                                                              y_vertex_cit->id());
          assert(connected);
          if (!connected || distance > 2) {
            continue;
          }
        }

        for (const auto& [ edge_label, 
                           edge_label_counter ] 
                         : graph_basic_statistics.edge_label_counter()) {
          // not add edge which does not exist in data graph
          if (!graph_basic_statistics.edge_type_counter().count(std::make_tuple(
                                    x_vertex_cit->label(),
                                    edge_label, 
                                    y_vertex_cit->label()))) {
            continue;
          }

          if (GUNDAM::HasEdge<const GraphPatternType>(
                  expand_tree_node.const_pattern().FindVertex(x_vertex_cit->id()), 
                  edge_label,
                  expand_tree_node.const_pattern().FindVertex(y_vertex_cit->id()))) {
            // this edge literal has duplicated with an existed edge
            continue;
          }
          // add to possible_edge_literals
          expanded_literal_info.emplace_back(
                              x_vertex_cit->id(),
                              y_vertex_cit->id(),
                              edge_label);
        }
      }
    }
  }

  if ( restriction.consider_constant_literal()
   &&(sorted_vertex_handle_set.empty() // has not specified vertex set 
   || sorted_vertex_handle_set.size() == 1)) { // has specified one vertexes
    // ########################
    // ##  constant literal  ##
    // ##  x.A = c           ##
    // ########################
    // enumerate x, A and c
    
    // enumerate x in the pattern
    for (auto x_vertex_cit = expand_tree_node.const_pattern().VertexBegin();
             !x_vertex_cit.IsDone(); 
              x_vertex_cit++){
      // enumerate x vertex
      PatternVertexConstHandleType x_vertex_handle = x_vertex_cit;
      if (!sorted_vertex_handle_set.empty()
       && !std::binary_search(sorted_vertex_handle_set.begin(),
                              sorted_vertex_handle_set.end(),
                                   x_vertex_handle)) {
        // specified vertex set but does not contain x_vertex_cit
        continue;
      }

      if (!std::binary_search(sort_new_vertex_id_set.begin(),
                              sort_new_vertex_id_set.end(),
                                     x_vertex_cit->id())) {
        // x is not contained in the new vertexes
        // move to the next vertex
        continue;
      }
      // iterator
      auto attr_it = graph_basic_statistics.legal_attr_set().find(x_vertex_cit->label());
      if (attr_it == graph_basic_statistics.legal_attr_set().cend()) {
        // this label does not have legal attr
        continue;
      }
      for (const auto& [attr_key, attr_value_container] : attr_it->second) {
        // enumerate A contained in x
        // A = attr_key
        for (const auto& c : attr_value_container) {
          // enumerate the c for x.A
          expanded_literal_info.emplace_back(x_vertex_cit->id(), 
                                             attr_key, 
                                             c.first.first,  // value_str
                                             c.first.second);// value_type
        }
      }
    }
  }

  // ######################################
  // ##  does not add ml literal in rhs  ##
  // ######################################
  if (restriction.consider_ml_literal()
   &&!rhs_literal_only // add ml literal only in lhs
   &&!ml_edge_label_set.empty()
   &&(sorted_vertex_handle_set.empty() // has not specified vertex set 
   || sorted_vertex_handle_set.size() == 2)) { // has specified one vertexes

    // ml literal
    // x_id_ - (module_url_: edge_label_) -> y_id_
    for (auto x_vertex_cit = expand_tree_node.const_pattern()
                                            .VertexBegin();
             !x_vertex_cit.IsDone(); 
              x_vertex_cit++) {
      // enumerate x vertex
      PatternVertexConstHandleType x_vertex_handle = x_vertex_cit;
      if (!sorted_vertex_handle_set.empty()
       && !std::binary_search(sorted_vertex_handle_set.begin(),
                              sorted_vertex_handle_set.end(),
                                   x_vertex_handle)) {
        // specified vertex set but does not contain x_vertex_cit
        continue;
      }
      for (auto y_vertex_cit = expand_tree_node.const_pattern()
                                              .VertexBegin();
               !y_vertex_cit.IsDone(); 
                y_vertex_cit++) {
        // not add self loop as ml literal
        if (x_vertex_cit->id() == y_vertex_cit->id()) 
          continue;
          
        PatternVertexConstHandleType y_vertex_handle = y_vertex_cit;
        if (!sorted_vertex_handle_set.empty()
         && !std::binary_search(sorted_vertex_handle_set.begin(),
                                sorted_vertex_handle_set.end(),
                                     y_vertex_handle)) {
          // specified vertex set but does not contain y_vertex_cit
          continue;
        }

        for (const auto& ml_edge_label : ml_edge_label_set) {
          // not add edge which does not exist in data graph
          if (!ml_edge_type_set.count(std::make_tuple(
                                    x_vertex_cit->label(),
                                    ml_edge_label, 
                                    y_vertex_cit->label()))) {
            continue;
          }

           auto  normal_edge_label_it  = ml_to_normal_edge_label.find(ml_edge_label);
          assert(normal_edge_label_it != ml_to_normal_edge_label.end());
          const auto& normal_edge_label = normal_edge_label_it->second;
          if (GUNDAM::HasEdge<const GraphPatternType>(
                  expand_tree_node.const_pattern().FindVertex(x_vertex_cit->id()), 
                  normal_edge_label,
                  expand_tree_node.const_pattern().FindVertex(y_vertex_cit->id()))) {
            assert(false); // ml literal should not duplicate with edge existed in pattern
            // this edge literal has duplicated with an existed edge
            continue;
          }
          // add to possible_edge_literals
          expanded_literal_info.emplace_back(
                              x_vertex_cit->id(),
                              y_vertex_cit->id(),
                              normal_edge_label,
                              "#"); // url
        }
      }
    }
  }
  
  if (rhs_literal_only
   && !restriction.specified_rhs_literal_set().empty()){
    // rhs has been specified
    gar::GraphAssociationRule<GraphPatternType,
                                 DataGraphType> 
         gar(expand_tree_node.const_pattern());
    for (auto literal_it  = expanded_literal_info.begin();
              literal_it != expanded_literal_info.end();){
      assert(gar.x_literal_set().Empty());
      // gar.x_literal_set().Clear();
      gar.ClearXLiteralSet();
      gar.AddX(*literal_it);
      assert(gar.x_literal_set().Count() == 1);
      auto stand_alone_info = (*gar.x_literal_set().begin())->stand_alone_info();
      bool is_specified = false;
      for (const auto& [specified_rhs_literal, support_set]
          : restriction.specified_rhs_literal_set()) {
        if (specified_rhs_literal == stand_alone_info){
          is_specified = true;
          break;
        }
      }
      if (!is_specified) {
        // is not contained in restriction.specified_rhs_literal_set()
        literal_it = expanded_literal_info.erase(literal_it);
        // gar.x_literal_set().Clear();
        gar.ClearXLiteralSet();
        continue;
      }
      // is contained in restriction.specified_rhs_literal_set()
      // gar.x_literal_set().Clear();
      gar.ClearXLiteralSet();
      literal_it++;
    }
  }
  return;
}

// assume that the vertex with the same label have the
// same kind of values
template <typename GraphPatternType,
          typename    DataGraphType>
void ExpandLiteral(
     ExpandTreeNode<GraphPatternType,
                       DataGraphType>& expand_tree_node,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
    const std::set<typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_edge_label_set, 
    const std::set<
        std::tuple<typename GUNDAM::VertexLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::VertexLabel<GraphPatternType>::type>>& ml_edge_type_set,
    const std::map<typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_to_normal_edge_label,
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
          std::vector<
           gar::LiteralInfo<GraphPatternType,
                               DataGraphType>>& expanded_literal_info) {

  std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type> vertex_handle_set;           

  _ExpandLiteral(expand_tree_node,
                vertex_handle_set,
           graph_basic_statistics,
                ml_edge_label_set, 
                 ml_edge_type_set,
          ml_to_normal_edge_label,
                      restriction, 
                            false,
            expanded_literal_info);

  return;
}

template <typename GraphPatternType,
          typename    DataGraphType>
void ExpandRhsLiteral(
    const ExpandTreeNode<GraphPatternType,
                            DataGraphType>& expand_tree_node,
    const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& vertex_handle_set,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
    const std::set<typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_edge_label_set, 
    const std::set<
        std::tuple<typename GUNDAM::VertexLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::VertexLabel<GraphPatternType>::type>>& ml_edge_type_set,
    const std::map<typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_to_normal_edge_label,
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
          std::vector<
           gar::LiteralInfo<GraphPatternType,
                               DataGraphType>>& expanded_literal_info) {

  _ExpandLiteral(expand_tree_node,
                vertex_handle_set,
           graph_basic_statistics,
                ml_edge_label_set, 
                 ml_edge_type_set,
          ml_to_normal_edge_label,
                      restriction, 
                             true,
            expanded_literal_info);

  return;
}

template <typename GraphPatternType,
          typename    DataGraphType>
void ExpandRhsLiteral(
    const ExpandTreeNode<GraphPatternType,
                            DataGraphType>& expand_tree_node,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
    const std::set<typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_edge_label_set, 
    const std::set<
        std::tuple<typename GUNDAM::VertexLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::VertexLabel<GraphPatternType>::type>>& ml_edge_type_set,
    const std::map<typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                   typename GUNDAM::  EdgeLabel<GraphPatternType>::type>& ml_to_normal_edge_label,
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
          std::vector<
           gar::LiteralInfo<GraphPatternType,
                               DataGraphType>>& expanded_literal_info) {

  std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type> vertex_handle_set;

  ExpandRhsLiteral(expand_tree_node,
                  vertex_handle_set,
             graph_basic_statistics,
                  ml_edge_label_set, 
                   ml_edge_type_set,
            ml_to_normal_edge_label,
                        restriction,
              expanded_literal_info);
  return;
}

template <typename LiteralType>
void ExpandableLhsLiterals(const std::vector<LiteralType>& all_possible_lhs_literals,
                           const std::vector<LiteralType>&          existed_literals,
                                 std::vector<LiteralType>&   expandable_lhs_literals) {
  // expandable_literals = all_possible_lhs_literals - existed_literals
  expandable_lhs_literals.clear();
  // the size of existed_literal is bounded,
  for (const auto& possible_literal : all_possible_lhs_literals) {
    bool existed_same = false;
    for (const auto& existed_literal : existed_literals) {
      if (possible_literal == existed_literal) {
        existed_same = true;
        break;
      }
    }
    if (!existed_same) {
      expandable_lhs_literals.emplace_back(possible_literal);
    }
  }
  assert((all_possible_lhs_literals.size()
                == existed_literals.size() 
         + expandable_lhs_literals.size())
      || (all_possible_lhs_literals.size() // the rhs literal contained in existed_literals might not be consiered in lhs
                == existed_literals.size() 
          + expandable_lhs_literals.size() - 1));
  return;
}

};  // namespace _gar_discover

}  // namespace grape

#endif  // _EXPAND_LITERAL_H