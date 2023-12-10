#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPAND_PATTERN_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_EXPAND_PATTERN_H_

#include "rule_discover/gar_discover/expand_tree.h"
#include "rule_discover/gar_discover/restriction.h"

#include "gundam/graph_type/small_graph.h"
#include "gundam/graph_type/graph.h"
#include "gundam/graph_type/large_graph.h"
#include "gundam/graph_type/large_graph2.h"

#include "gundam/type_getter/vertex_handle.h"
#include "gundam/type_getter/vertex_label.h"
#include "gundam/type_getter/vertex_id.h"
#include "gundam/type_getter/edge_handle.h"
#include "gundam/type_getter/edge_label.h"
#include "gundam/type_getter/edge_id.h"

#include "gundam/tool/isolate_vertex.h"
#include "gundam/tool/topological/path/is_path.h"
#include "gundam/tool/topological/star/is_star.h"
#include "gundam/tool/topological/star/decompose_star_to_path.h"
#include "gundam/tool/diameter.h"
#include "gundam/tool/max_id.h"
#include "gundam/tool/connected_component.h"

#include "gundam/graph_statistics/graph_basic_statistics.h"

namespace grape{

namespace _gar_discover {

template <typename GraphPatternType, 
          typename    DataGraphType>
inline void AddNewEdgeInPattern(
            const ExpandTreeNode<GraphPatternType,
                                    DataGraphType>& expand_node,
 const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
  const typename GUNDAM::EdgeID<GraphPatternType>::type new_edge_id,
                 ExpandTreeLevel<GraphPatternType,
                                    DataGraphType> &expand_level) {

  auto& expand_node_pattern = expand_node.pattern();
                              
  auto pattern_size = expand_node_pattern.CountVertex();

  bool exist_duplicated_edge = false;

  using VertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using VertexIDType     = typename GUNDAM::    VertexID<GraphPatternType>::type;
  using EdgeLabelType    = typename GUNDAM::   EdgeLabel<GraphPatternType>::type;

  for (auto src_vertex_iter = expand_node_pattern.VertexBegin();
           !src_vertex_iter.IsDone();
            src_vertex_iter++) {
    for (auto dst_vertex_iter = expand_node_pattern.VertexBegin();
             !dst_vertex_iter.IsDone(); 
              dst_vertex_iter++) {
      // not add self loop
      if (src_vertex_iter->id() 
       == dst_vertex_iter->id()) 
        continue;
      if (!exist_duplicated_edge) {
        // to find whether src_vertex_iter and dst_vertex_iter
        // are connected
        bool connected = false;
        for (auto out_edge_it = src_vertex_iter->OutEdgeBegin();
                 !out_edge_it.IsDone();
                  out_edge_it++)  {
          if (out_edge_it->dst_handle() != dst_vertex_iter) {
            continue;
          }
          connected = true;
          break;
        }
        if (connected) {
          continue;
        }
        for (auto in_edge_it = src_vertex_iter->InEdgeBegin();
                 !in_edge_it.IsDone();
                  in_edge_it++)  {
          if (in_edge_it->src_handle() != dst_vertex_iter) {
            continue;
          }
          connected = true;
          break;
        }
        if (connected) {
          continue;
        }
      }
      for (const auto& [edge_label, edge_label_counter] 
           : graph_basic_statistics.edge_label_counter()) {
        // not add edge which does not exist in data graph
        if (!graph_basic_statistics.edge_type_counter().count(std::tuple(
                                  src_vertex_iter->label(),
                                  edge_label, 
                                  dst_vertex_iter->label()))) {
          continue;
        }
        GraphPatternType expand_pattern(expand_node_pattern);
        expand_pattern.AddEdge(src_vertex_iter->id(), 
                               dst_vertex_iter->id(),
                               edge_label, 
                               new_edge_id);
        expand_level.AddExpandTreeNode(
                            expand_level.size(),
                            expand_pattern,
                            // does not have new vertexes
                            expand_node.const_rhs_literals(),
                            expand_node.const_lhs_literals());
      }
    }
  }
  return;
}

template <typename GraphPatternType, 
          typename    DataGraphType>
inline void AddNewEdgeInPattern(GraphPatternType &graph_pattern,
 const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
  const typename GUNDAM::EdgeID<GraphPatternType>::type new_edge_id,
                 std::vector<GraphPatternType>& expand_pattern_vec) {

  auto& expand_node_pattern = graph_pattern;
                              
  auto pattern_size = expand_node_pattern.CountVertex();

  bool exist_duplicated_edge = false;

  using VertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using VertexIDType     = typename GUNDAM::    VertexID<GraphPatternType>::type;
  using EdgeLabelType    = typename GUNDAM::   EdgeLabel<GraphPatternType>::type;

  for (auto src_vertex_iter = expand_node_pattern.VertexBegin();
           !src_vertex_iter.IsDone();
            src_vertex_iter++) {
    for (auto dst_vertex_iter = expand_node_pattern.VertexBegin();
             !dst_vertex_iter.IsDone(); 
              dst_vertex_iter++) {
      // not add self loop
      if (src_vertex_iter->id() 
       == dst_vertex_iter->id()) 
        continue;
      if (!exist_duplicated_edge) {
        // to find whether src_vertex_iter and dst_vertex_iter
        // are connected
        bool connected = false;
        for (auto out_edge_it = src_vertex_iter->OutEdgeBegin();
                 !out_edge_it.IsDone();
                  out_edge_it++)  {
          if (out_edge_it->dst_handle() != dst_vertex_iter) {
            continue;
          }
          connected = true;
          break;
        }
        if (connected) {
          continue;
        }
        for (auto in_edge_it = src_vertex_iter->InEdgeBegin();
                 !in_edge_it.IsDone();
                  in_edge_it++)  {
          if (in_edge_it->src_handle() != dst_vertex_iter) {
            continue;
          }
          connected = true;
          break;
        }
        if (connected) {
          continue;
        }
      }
      for (const auto& [edge_label, edge_label_counter] 
           : graph_basic_statistics.edge_label_counter()) {
        // not add edge which does not exist in data graph
        if (!graph_basic_statistics.edge_type_counter().count(std::tuple(
                                  src_vertex_iter->label(),
                                  edge_label, 
                                  dst_vertex_iter->label()))) {
          continue;
        }
        GraphPatternType expand_pattern(expand_node_pattern);
        expand_pattern.AddEdge(src_vertex_iter->id(), 
                               dst_vertex_iter->id(),
                               edge_label, 
                               new_edge_id);
        expand_pattern_vec.emplace_back(expand_pattern);
      }
    }
  }
  return;
}




template <typename GraphPatternType, 
          typename    DataGraphType>
inline void AddNewEdgeOutPattern(
    const ExpandTreeNode<GraphPatternType,
                            DataGraphType> &expand_node,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
      const typename GUNDAM::EdgeID<GraphPatternType>::type new_edge_id, 
                                      bool only_expand_link,
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
    ExpandTreeLevel<GraphPatternType,
                       DataGraphType> &expand_level) {

  using     VertexIDType = typename GUNDAM::    VertexID<GraphPatternType>::type;
  using    EdgeLabelType = typename GUNDAM::   EdgeLabel<GraphPatternType>::type;
  using  VertexLabelType = typename GUNDAM:: VertexLabel<GraphPatternType>::type;
  using VertexHandleType = typename GUNDAM::VertexHandle<const GraphPatternType>::type;

  using EdgeTypeType = std::tuple<VertexLabelType,
                                    EdgeLabelType,
                                  VertexLabelType>;

  auto& expand_node_pattern = expand_node.pattern();

  const VertexIDType kNewVertexId = GUNDAM::MaxVertexId(expand_node_pattern) + 1;

  assert(!only_expand_link
       || GUNDAM::IsPath<true>(expand_node_pattern));

  // set of vertexes that can expand outward edges
  std::set<VertexHandleType> handle_set_can_expand_from;

  if (only_expand_link || restriction.pattern_is_link()
                       || restriction.pattern_is_star()) {
    // has specified pattern topologic, cannot expand outward/inward
    // edge from arbitiary vertex
    std::vector<GraphPatternType> 
      connected_components(GUNDAM::ConnectedComponent(expand_node_pattern));
    assert((!restriction.gcr() && connected_components.size() == 1)
         ||( restriction.gcr() && connected_components.size() == 2));
    for (const auto& connected_component 
                   : connected_components) {
      if (only_expand_link || restriction.pattern_is_link()) {
        auto [src_handle, 
              dst_handle] = GUNDAM::PathEndPoints<true>(connected_component);
        if (!src_handle) {
          assert(!dst_handle);
          // the pattern is not link before expanded
          return;
        }
        assert((src_handle != dst_handle)
             || connected_component.CountVertex() == 1);
        handle_set_can_expand_from.emplace(expand_node_pattern.FindVertex(src_handle->id()));
        handle_set_can_expand_from.emplace(expand_node_pattern.FindVertex(dst_handle->id()));
        assert(!handle_set_can_expand_from.empty());
        continue;
      }
      assert(restriction.pattern_is_star());
      const auto [end_vertex_handle_set,
              central_vertex_handle] = restriction.bidirectional_path()?
                                       GUNDAM::StarEndPoints< true>(connected_component)
                                     : GUNDAM::StarEndPoints<false>(connected_component);
      assert(end_vertex_handle_set.size() >= 2);
      assert(end_vertex_handle_set.size() == 2 || central_vertex_handle);
      if (!central_vertex_handle) {
        // is link, all vertexes can be legal central vertex if consider biditional path
        assert(( restriction.bidirectional_path() && GUNDAM::IsPath< true>(connected_component))
            || (!restriction.bidirectional_path() && GUNDAM::IsPath<false>(connected_component)));
        if (!restriction.specified_path_num_limit()
          ||(restriction.          path_num_limit() > 2)) {
          // has not reached the path num limit, add them all
          for (auto vertex_it = connected_component.VertexBegin();
                   !vertex_it.IsDone();
                    vertex_it++) {
            assert(expand_node_pattern
                              .FindVertex(vertex_it->id()));
            auto [ ends_handle_set_it,
                   ends_handle_set_ret ] 
                      = handle_set_can_expand_from.emplace(expand_node_pattern
                                                          .FindVertex(vertex_it->id()));
            assert(ends_handle_set_ret);
          }
          continue;
        }
        assert(restriction.path_num_limit() == 1
            || restriction.path_num_limit() == 2);
        // has reached the path num limit
        // add only the two end vertexes
        for (const auto& end_vertex_handle 
                       : end_vertex_handle_set) {
          assert(expand_node_pattern
                .FindVertex(end_vertex_handle->id()));
          handle_set_can_expand_from.emplace(expand_node_pattern
                                            .FindVertex(end_vertex_handle->id()));
        }
        continue;
      }
      // add all end vertexes and the central vertex
      for (const auto& vertex_handle : end_vertex_handle_set)  {
        assert(expand_node_pattern
              .FindVertex(vertex_handle->id()));
        if ( vertex_handle == central_vertex_handle ) {
          assert(false);
          continue;
        }
        handle_set_can_expand_from.emplace(expand_node_pattern
                                          .FindVertex(vertex_handle->id()));
      }
      if (restriction.specified_path_num_limit()
       &&(restriction.          path_num_limit() <= end_vertex_handle_set.size())) {
        // has reached the path num limit, does not added central vertex
        continue;
      }
      // has not reached the path num limit, add the central vertex
      assert(expand_node_pattern
            .FindVertex(central_vertex_handle->id()));
      auto [ ends_handle_set_it,
             ends_handle_set_ret ] 
           = handle_set_can_expand_from.emplace(expand_node_pattern
                                               .FindVertex(central_vertex_handle->id()));
      assert(ends_handle_set_ret);
    }
    assert(handle_set_can_expand_from.size() >= 2);
  }

  std::string expand_node_pattern_str;
  expand_node_pattern_str << expand_node_pattern;
  util::Debug("expand_node_pattern_str: " + expand_node_pattern_str);
  util::Debug("handle_set_can_expand_from.size(): " + std::to_string(handle_set_can_expand_from.size()));
  for (const auto& end_handle : handle_set_can_expand_from) {
    util::Debug("\tend_handle->id(): " + std::to_string(end_handle->id()));
  }

  #ifndef NDEBUG
  for (const auto& vertex_handle 
                 : handle_set_can_expand_from) {
    assert(expand_node_pattern
          .FindVertex(vertex_handle->id()) == vertex_handle);
  }
  #endif // NDEBUG

  for (auto vertex_iter = expand_node_pattern.VertexBegin(); 
           !vertex_iter.IsDone();
            vertex_iter++) {
    if (!handle_set_can_expand_from.empty()) {
      if (handle_set_can_expand_from.find(vertex_iter) 
       == handle_set_can_expand_from.end()) {
        // cannot expand edge from this vertex
        continue;
      }
    }
    util::Debug("Expand from vertex id: " + std::to_string(vertex_iter->id()));
    for (const auto& [vertex_label, vertex_label_counter] 
           : graph_basic_statistics.vertex_label_counter()) {
      for (const auto& [edge_label,   edge_label_counter]    
             : graph_basic_statistics.edge_label_counter()) {
        // not add edge which does not exist in data graph
        auto expand_pattern_func = [&](const auto& is_output) -> void {
          if (is_output) {
            util::Debug("Expand output edge");
          }
          else {
            util::Debug("Expand input edge");
          }
          const auto kNewEdgeType 
                   = is_output ? std::tuple(vertex_iter->label(),
                                              edge_label, 
                                            vertex_label)
                               : std::tuple(vertex_label,
                                              edge_label, 
                                            vertex_iter->label());
          if (!graph_basic_statistics.edge_type_counter()
                                     .count(kNewEdgeType)) {
            return;
          }
          // and an out edge to new vertex 
          for (auto out_edge_it = vertex_iter->OutEdgeBegin();
                   !out_edge_it.IsDone();
                    out_edge_it++) {
            EdgeTypeType existed_edge_type(out_edge_it->src_handle()->label(),
                                           out_edge_it->label(), 
                                           out_edge_it->dst_handle()->label());
            if (graph_basic_statistics.EdgeTypeConnectedTo(existed_edge_type, false,
                                                                kNewEdgeType, !is_output)) {
              continue;
            } 
            // does not connected to existed edge types
            return;
          }
          for (auto in_edge_it = vertex_iter->InEdgeBegin();
                   !in_edge_it.IsDone();
                    in_edge_it++) {
            EdgeTypeType existed_edge_type(in_edge_it->src_handle()->label(),
                                           in_edge_it->label(), 
                                           in_edge_it->dst_handle()->label());
            if (graph_basic_statistics.EdgeTypeConnectedTo(existed_edge_type, true,
                                                                kNewEdgeType, !is_output)) {
              continue;
            }
            // does not connected to existed edge types
            return;
          }
          GraphPatternType expand_pattern(expand_node_pattern);
          expand_pattern.AddVertex(kNewVertexId,
                                       vertex_label);
          if (is_output) {
            expand_pattern.AddEdge(vertex_iter->id(), 
                               kNewVertexId,
                                     edge_label, 
                                 new_edge_id);
          }
          else {
            expand_pattern.AddEdge(kNewVertexId, 
                                       vertex_iter->id(),
                                         edge_label, 
                                     new_edge_id);
          }
          std::string expand_pattern_str;
          expand_pattern_str << expand_pattern;
          util::Debug("expand_pattern_str: " 
                     + expand_pattern_str);
          if (restriction.has_diameter_limit()
           &&(restriction.diameter_limit() < GUNDAM::Diameter<true>(expand_pattern))) {
            // reached the diameter_limit, would not expand
            util::Debug("not satisfy diameter_limit");
            return;
          }
          if (!restriction.bidirectional_path()) {
            // to find whether this is a star with outward edges

            std::vector<GraphPatternType> 
                       connected_components(GUNDAM::ConnectedComponent(expand_pattern));
            for (const auto& connected_component : connected_components) {
              const auto [end_vertex_handle_set,
                      central_vertex_handle] = GUNDAM::StarEndPoints<false>(connected_component);
              if (end_vertex_handle_set.empty()) {
                // is not star
                util::Debug("is not legal pattern when BidirectionalPath set to false");
                return;
              }
              if (!central_vertex_handle) {
                // is a legal single directed path
                continue;
              }
              // find the central_vertex_handle, find whether the 
              // edges are all outward
              if (central_vertex_handle->CountInEdge() > 0) {
                // has input edge, illegal
                util::Debug("is not legal pattern when BidirectionalPath set to false");
                return;
              }
            }
          }
          if (restriction.specified_path_length_limit()) {
            std::vector<GraphPatternType> 
                   connected_components(GUNDAM::ConnectedComponent(expand_pattern));
            assert(connected_components.size() == 2);
            for (auto& connected_component : connected_components) {
              if (restriction.bidirectional_path()) {
                if (GUNDAM::IsPath<true>(connected_component)) {
                  if (connected_component.CountEdge() <= 2 * restriction.path_length_limit()) {
                    // each vertex can be the central vertex
                    // satisfy the path_length_limit
                    continue;
                  }
                  // does not satisfy the path_length_limit
                  // for regarding each vertex as the central vertex
                  util::Debug("reached specified_path_length_limit");
                  return;
                }
              }
              auto path_set = restriction.bidirectional_path()?
                              GUNDAM::DecomposeStarToPath< true>(connected_component)
                            : GUNDAM::DecomposeStarToPath<false>(connected_component);
              std::string connected_component_str;
              connected_component_str << connected_component;
              util::Debug("connected_component_str: " + connected_component_str);
              assert(path_set.size() <= 2);
              for (const auto& path : path_set) {
                std::string path_str;
                path_str << path;
                util::Debug("\tpath_str: " + path_str);
                if (path.CountEdge() <= restriction.path_length_limit()) {
                  // satisfy the path_length_limit
                  continue;
                }
                // does not satisfy the path_length_limit
                util::Debug("reached specified_path_length_limit");
                return;
              }
            }
          }
          // legal, can be added in expand_level
          expand_level.AddExpandTreeNode(
                          expand_level.size(),
                          expand_pattern,
                          kNewVertexId,
                          expand_node.const_rhs_literals(),
                          expand_node.const_lhs_literals());
        };

        expand_pattern_func( true); // output
        expand_pattern_func(false); //  input
      }
    }
  }
  return;
}

template <typename GraphPatternType,
          typename    DataGraphType>
inline void ExpandPattern(
               const ExpandTreeNode<GraphPatternType,
                                       DataGraphType> &expand_node,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
    const
    typename GUNDAM::EdgeID<GraphPatternType>::type new_edge_id, 
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
                                       bool only_expand_link,
            ExpandTreeLevel<GraphPatternType,
                               DataGraphType> &expand_level) {
  util::Debug("only_expand_link: " + std::to_string(only_expand_link));
  util::Debug("restriction.pattern_is_tree(): " + std::to_string(restriction.pattern_is_tree()));
  util::Debug("restriction.pattern_is_star(): " + std::to_string(restriction.pattern_is_star()));
  util::Debug("restriction.pattern_is_link(): " + std::to_string(restriction.pattern_is_link()));
  if (!only_expand_link
   && !restriction.pattern_is_tree()
   && !restriction.pattern_is_star()
   && !restriction.pattern_is_link()
   && expand_node.const_rhs_literals().size() != 0){
    // AddNewEdgeInPattern won't increase possible
    // literals, if does not have possible literal
    // before, should not added in edge
    AddNewEdgeInPattern(expand_node,
                        graph_basic_statistics, 
                        new_edge_id, 
                        expand_level);
  }
  util::Debug("here!");
  if (restriction.has_pattern_vertex_limit()
   &&(restriction.    pattern_vertex_limit() 
    < expand_node.    pattern().CountVertex())) {
    return;
  }
  AddNewEdgeOutPattern(expand_node,
                       graph_basic_statistics, 
                    new_edge_id,
                    only_expand_link,
                    restriction,
                    expand_level);
  return;
}




template <typename GraphPatternType, 
          typename    DataGraphType>
inline void AddNewEdgeOutPattern(GraphPatternType &graph_pattern,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
      const typename GUNDAM::EdgeID<GraphPatternType>::type new_edge_id, 
                                      bool only_expand_link,
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
    std::vector<GraphPatternType>& expand_pattern_vec) {

  using     VertexIDType = typename GUNDAM::    VertexID<GraphPatternType>::type;
  using    EdgeLabelType = typename GUNDAM::   EdgeLabel<GraphPatternType>::type;
  using  VertexLabelType = typename GUNDAM:: VertexLabel<GraphPatternType>::type;
  using VertexHandleType = typename GUNDAM::VertexHandle<const GraphPatternType>::type;

  using EdgeTypeType = std::tuple<VertexLabelType,
                                    EdgeLabelType,
                                  VertexLabelType>;

  const auto& expand_node_pattern = graph_pattern;

  const VertexIDType kNewVertexId = GUNDAM::MaxVertexId(expand_node_pattern) + 1;

  assert(!only_expand_link
       || GUNDAM::IsPath<true>(expand_node_pattern));

  // set of vertexes that can expand outward edges
  std::set<VertexHandleType> handle_set_can_expand_from;

  if (only_expand_link || restriction.pattern_is_link()
                       || restriction.pattern_is_star()) {
    // has specified pattern topologic, cannot expand outward/inward
    // edge from arbitiary vertex
    std::vector<GraphPatternType> 
      connected_components(GUNDAM::ConnectedComponent(expand_node_pattern));
    assert((!restriction.gcr() && connected_components.size() == 1)
         ||( restriction.gcr() && connected_components.size() == 2));
    for (const auto& connected_component 
                   : connected_components) {
      if (only_expand_link || restriction.pattern_is_link()) {
        auto [src_handle, 
              dst_handle] = GUNDAM::PathEndPoints<true>(connected_component);
        if (!src_handle) {
          assert(!dst_handle);
          // the pattern is not link before expanded
          return;
        }
        assert((src_handle != dst_handle)
             || connected_component.CountVertex() == 1);
        handle_set_can_expand_from.emplace(expand_node_pattern.FindVertex(src_handle->id()));
        handle_set_can_expand_from.emplace(expand_node_pattern.FindVertex(dst_handle->id()));
        assert(!handle_set_can_expand_from.empty());
        continue;
      }
      assert(restriction.pattern_is_star());
      const auto [end_vertex_handle_set,
              central_vertex_handle] = restriction.bidirectional_path()?
                                       GUNDAM::StarEndPoints< true>(connected_component)
                                     : GUNDAM::StarEndPoints<false>(connected_component);
      assert(end_vertex_handle_set.size() >= 2);
      assert(end_vertex_handle_set.size() == 2 || central_vertex_handle);
      if (!central_vertex_handle) {
        // is link, all vertexes can be legal central vertex if consider biditional path
        assert(( restriction.bidirectional_path() && GUNDAM::IsPath< true>(connected_component))
            || (!restriction.bidirectional_path() && GUNDAM::IsPath<false>(connected_component)));
        if (!restriction.specified_path_num_limit()
          ||(restriction.          path_num_limit() > 2)) {
          // has not reached the path num limit, add them all
          for (auto vertex_it = connected_component.VertexBegin();
                   !vertex_it.IsDone();
                    vertex_it++) {
            assert(expand_node_pattern
                              .FindVertex(vertex_it->id()));
            auto [ ends_handle_set_it,
                   ends_handle_set_ret ] 
                      = handle_set_can_expand_from.emplace(expand_node_pattern
                                                          .FindVertex(vertex_it->id()));
            assert(ends_handle_set_ret);
          }
          continue;
        }
        assert(restriction.path_num_limit() == 1
            || restriction.path_num_limit() == 2);
        // has reached the path num limit
        // add only the two end vertexes
        for (const auto& end_vertex_handle 
                       : end_vertex_handle_set) {
          assert(expand_node_pattern
                .FindVertex(end_vertex_handle->id()));
          handle_set_can_expand_from.emplace(expand_node_pattern
                                            .FindVertex(end_vertex_handle->id()));
        }
        continue;
      }
      // add all end vertexes and the central vertex
      for (const auto& vertex_handle : end_vertex_handle_set)  {
        assert(expand_node_pattern
              .FindVertex(vertex_handle->id()));
        if ( vertex_handle == central_vertex_handle ) {
          assert(false);
          continue;
        }
        handle_set_can_expand_from.emplace(expand_node_pattern
                                          .FindVertex(vertex_handle->id()));
      }
      if (restriction.specified_path_num_limit()
       &&(restriction.          path_num_limit() <= end_vertex_handle_set.size())) {
        // has reached the path num limit, does not added central vertex
        continue;
      }
      // has not reached the path num limit, add the central vertex
      assert(expand_node_pattern
            .FindVertex(central_vertex_handle->id()));
      auto [ ends_handle_set_it,
             ends_handle_set_ret ] 
           = handle_set_can_expand_from.emplace(expand_node_pattern
                                               .FindVertex(central_vertex_handle->id()));
      assert(ends_handle_set_ret);
    }
    assert(handle_set_can_expand_from.size() >= 2);
  }

  std::string expand_node_pattern_str;
  expand_node_pattern_str << expand_node_pattern;
  util::Debug("expand_node_pattern_str: " + expand_node_pattern_str);
  util::Debug("handle_set_can_expand_from.size(): " + std::to_string(handle_set_can_expand_from.size()));
  for (const auto& end_handle : handle_set_can_expand_from) {
    util::Debug("\tend_handle->id(): " + std::to_string(end_handle->id()));
  }

  #ifndef NDEBUG
  for (const auto& vertex_handle 
                 : handle_set_can_expand_from) {
    assert(expand_node_pattern
          .FindVertex(vertex_handle->id()) == vertex_handle);
  }
  #endif // NDEBUG

  for (auto vertex_iter = expand_node_pattern.VertexBegin(); 
           !vertex_iter.IsDone();
            vertex_iter++) {
    if (!handle_set_can_expand_from.empty()) {
      if (handle_set_can_expand_from.find(vertex_iter) 
       == handle_set_can_expand_from.end()) {
        // cannot expand edge from this vertex
        continue;
      }
    }
    util::Debug("Expand from vertex id: " + std::to_string(vertex_iter->id()));
    for (const auto& [vertex_label, vertex_label_counter] 
           : graph_basic_statistics.vertex_label_counter()) {
      for (const auto& [edge_label,   edge_label_counter]    
             : graph_basic_statistics.edge_label_counter()) {
        // not add edge which does not exist in data graph
        auto expand_pattern_func = [&](const auto& is_output) -> void {
          if (is_output) {
            util::Debug("Expand output edge");
          }
          else {
            util::Debug("Expand input edge");
          }
          const auto kNewEdgeType 
                   = is_output ? std::tuple(vertex_iter->label(),
                                              edge_label, 
                                            vertex_label)
                               : std::tuple(vertex_label,
                                              edge_label, 
                                            vertex_iter->label());
          if (!graph_basic_statistics.edge_type_counter()
                                     .count(kNewEdgeType)) {
            return;
          }
          // and an out edge to new vertex 
          for (auto out_edge_it = vertex_iter->OutEdgeBegin();
                   !out_edge_it.IsDone();
                    out_edge_it++) {
            EdgeTypeType existed_edge_type(out_edge_it->src_handle()->label(),
                                           out_edge_it->label(), 
                                           out_edge_it->dst_handle()->label());
            if (graph_basic_statistics.EdgeTypeConnectedTo(existed_edge_type, false,
                                                                kNewEdgeType, !is_output)) {
              continue;
            } 
            // does not connected to existed edge types
            return;
          }
          for (auto in_edge_it = vertex_iter->InEdgeBegin();
                   !in_edge_it.IsDone();
                    in_edge_it++) {
            EdgeTypeType existed_edge_type(in_edge_it->src_handle()->label(),
                                           in_edge_it->label(), 
                                           in_edge_it->dst_handle()->label());
            if (graph_basic_statistics.EdgeTypeConnectedTo(existed_edge_type, true,
                                                                kNewEdgeType, !is_output)) {
              continue;
            }
            // does not connected to existed edge types
            return;
          }
          GraphPatternType expand_pattern(expand_node_pattern);
          expand_pattern.AddVertex(kNewVertexId,
                                       vertex_label);
          if (is_output) {
            expand_pattern.AddEdge(vertex_iter->id(), 
                               kNewVertexId,
                                     edge_label, 
                                 new_edge_id);
          }
          else {
            expand_pattern.AddEdge(kNewVertexId, 
                                       vertex_iter->id(),
                                         edge_label, 
                                     new_edge_id);
          }
          std::string expand_pattern_str;
          expand_pattern_str << expand_pattern;
          util::Debug("expand_pattern_str: " 
                     + expand_pattern_str);
          if (restriction.has_diameter_limit()
           &&(restriction.diameter_limit() < GUNDAM::Diameter<true>(expand_pattern))) {
            // reached the diameter_limit, would not expand
            util::Debug("not satisfy diameter_limit");
            return;
          }
          if (!restriction.bidirectional_path()) {
            // to find whether this is a star with outward edges

            std::vector<GraphPatternType> 
                       connected_components(GUNDAM::ConnectedComponent(expand_pattern));
            for (const auto& connected_component : connected_components) {
              const auto [end_vertex_handle_set,
                      central_vertex_handle] = GUNDAM::StarEndPoints<false>(connected_component);
              if (end_vertex_handle_set.empty()) {
                // is not star
                util::Debug("is not legal pattern when BidirectionalPath set to false");
                return;
              }
              if (!central_vertex_handle) {
                // is a legal single directed path
                continue;
              }
              // find the central_vertex_handle, find whether the 
              // edges are all outward
              if (central_vertex_handle->CountInEdge() > 0) {
                // has input edge, illegal
                util::Debug("is not legal pattern when BidirectionalPath set to false");
                return;
              }
            }
          }
          if (restriction.specified_path_length_limit()) {
            std::vector<GraphPatternType> 
                   connected_components(GUNDAM::ConnectedComponent(expand_pattern));
            assert(connected_components.size() == 2);
            for (auto& connected_component : connected_components) {
              if (restriction.bidirectional_path()) {
                if (GUNDAM::IsPath<true>(connected_component)) {
                  if (connected_component.CountEdge() <= 2 * restriction.path_length_limit()) {
                    // each vertex can be the central vertex
                    // satisfy the path_length_limit
                    continue;
                  }
                  // does not satisfy the path_length_limit
                  // for regarding each vertex as the central vertex
                  util::Debug("reached specified_path_length_limit");
                  return;
                }
              }
              auto path_set = restriction.bidirectional_path()?
                              GUNDAM::DecomposeStarToPath< true>(connected_component)
                            : GUNDAM::DecomposeStarToPath<false>(connected_component);
              std::string connected_component_str;
              connected_component_str << connected_component;
              util::Debug("connected_component_str: " + connected_component_str);
              assert(path_set.size() <= 2);
              for (const auto& path : path_set) {
                std::string path_str;
                path_str << path;
                util::Debug("\tpath_str: " + path_str);
                if (path.CountEdge() <= restriction.path_length_limit()) {
                  // satisfy the path_length_limit
                  continue;
                }
                // does not satisfy the path_length_limit
                util::Debug("reached specified_path_length_limit");
                return;
              }
            }
          }
          // legal, can be added in expand_level
          expand_pattern_vec.emplace_back(expand_pattern);
        };

        expand_pattern_func( true); // output
        expand_pattern_func(false); //  input
      }
    }
  }
  return;
}




template <typename GraphPatternType,
          typename    DataGraphType>
inline std::vector<GraphPatternType> ExpandPattern(GraphPatternType &graph_pattern,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
    const
    typename GUNDAM::EdgeID<GraphPatternType>::type new_edge_id,
          const Restriction<GraphPatternType,
                               DataGraphType>& restriction,
                                       bool only_expand_link) {
  util::Debug("only_expand_link: " + std::to_string(only_expand_link));
  util::Debug("restriction.pattern_is_tree(): " + std::to_string(restriction.pattern_is_tree()));
  util::Debug("restriction.pattern_is_star(): " + std::to_string(restriction.pattern_is_star()));
  util::Debug("restriction.pattern_is_link(): " + std::to_string(restriction.pattern_is_link()));
  std::vector<GraphPatternType> expand_pattern_vec;
  if (!only_expand_link
   && !restriction.pattern_is_tree()
   && !restriction.pattern_is_star()
   && !restriction.pattern_is_link()){
    AddNewEdgeInPattern(graph_pattern,
                        graph_basic_statistics, 
                        new_edge_id, 
                        expand_pattern_vec);
  }
  util::Debug("here!");
  if (restriction.has_pattern_vertex_limit()
   &&(restriction.    pattern_vertex_limit() 
    < graph_pattern.           CountVertex())) {
    return expand_pattern_vec;
  }
  AddNewEdgeOutPattern(graph_pattern,
                       graph_basic_statistics, 
                    new_edge_id,
                    only_expand_link,
                    restriction,
                    expand_pattern_vec);
  return expand_pattern_vec;
}




} // namespace _gar_discover

} // namespace grape

#endif // _EXPAND_PATTERN_H
