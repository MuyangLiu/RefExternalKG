#ifndef _UNIQUE_PATTERN_H
#define _UNIQUE_PATTERN_H

#include <omp.h>

#include "gundam/graph_type/small_graph.h"
#include "gundam/graph_type/graph.h"
#include "gundam/graph_type/large_graph.h"
#include "gundam/graph_type/large_graph2.h"

#include "gundam/algorithm/vf2.h"

#include "gundam/tool/same_pattern.h"
#include "gundam/tool/topological/path/is_path.h"
#include "gundam/tool/has_edge.h"
#include "gundam/tool/deduplicate_patterns/dfs_code.h"
#include "gundam/tool/deduplicate_patterns/deduplicate_patterns.h"

#include "gundam/component/disjoint_set.h"

#include "gundam/type_getter/vertex_handle.h"
#include "gundam/type_getter/vertex_label.h"
#include "gundam/type_getter/edge_label.h"

#include "rule_discover/gar_discover/gar_dfs_code.h"
#include "rule_discover/gar_discover/expand_literal.h"
#include "rule_discover/gar_discover/restriction.h"

#include "gar/same_gar.h"

namespace grape{

namespace _gar_discover {

template<typename GraphPatternType,
         typename    DataGraphType>
inline void MergeExpandedPatterns(
                      ExpandTreeLevel<GraphPatternType,
                                         DataGraphType>& expand_level,
                      ExpandTreeLevel<GraphPatternType,
                                         DataGraphType>& before_expand_level,
                                    std::vector<size_t>& expand_from_pattern,
              const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_basic_statistics,
  const std::set<typename GUNDAM::  EdgeLabel<GraphPatternType>::type>&  ml_edge_label_set, 
  const std::set<
      std::tuple<typename GUNDAM::VertexLabel<GraphPatternType>::type,
                 typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                 typename GUNDAM::VertexLabel<GraphPatternType>::type>>& ml_edge_type_set,
  const std::map<typename GUNDAM::  EdgeLabel<GraphPatternType>::type,
                 typename GUNDAM::  EdgeLabel<GraphPatternType>::type>&  ml_to_normal_edge_label,
                            const Restriction<GraphPatternType,
                                                 DataGraphType>& restriction) {

  using ExpandTreeLevelType 
      = ExpandTreeLevel<GraphPatternType,
                           DataGraphType>;

  using DfsCodeType = GUNDAM::DfsCode<GraphPatternType>;

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using GarType = gar::GraphAssociationRule<GraphPatternType, 
                                               DataGraphType>;

  using LiteralInfoType = gar::LiteralInfo<GraphPatternType, 
                                              DataGraphType>;

  // patterns with same dfs code are same
  std::map<DfsCodeType, std::set<size_t>> dfs_code_map;

  // ########################################################
  // ##  first generate DfsCode set for each expanded      ##
  // ##  pattern that has either non-empty new_vertex_set  ##
  // ##  or non-empty legal literal set                    ##
  // ########################################################
  util::Info("begin to par generate dfs code");
  omp_lock_t     generate_dfs_code_lock; 
  omp_init_lock(&generate_dfs_code_lock);

  #pragma omp parallel for schedule(dynamic) 
  for (size_t expand_node_id = 0; 
              expand_node_id < expand_level.size();
              expand_node_id++) {
    auto& expand_node = expand_level.node(expand_node_id);
    if (expand_node.new_vertexes_id().empty()
     && expand_node.   rhs_literals().empty()) {
      // this pattern does not have legal rhs, and also
      // did not add new vertex, would be used to prune
      // expanded pattern later
      continue;
    }
    std::vector<DfsCodeType> dfs_code_set_for_pattern 
             = GUNDAM::GetDFSCode(expand_node.pattern());
    assert(!dfs_code_set_for_pattern.empty());

    omp_set_lock(&generate_dfs_code_lock);
    for (const auto& dfs_code : dfs_code_set_for_pattern) {
      dfs_code_map[dfs_code].emplace(expand_node_id);
    }
    omp_unset_lock(&generate_dfs_code_lock);
  }
  util::Info("dfs code generated");

  // ###################################################
  // ##  then generate DfsCode set for each expanded  ##
  // ##  pattern that has both empty new_vertex_set   ##
  // ##  and empty legal literal set for pruning      ##
  // ###################################################
  util::Info("begin to generate dfs code for pruning");
  for (size_t expand_node_id = 0; 
              expand_node_id < expand_level.size();
              expand_node_id++) {
    auto& expand_node = expand_level.node(expand_node_id);
    if (!expand_node.new_vertexes_id().empty()
     || !expand_node.   rhs_literals().empty()) {
      // this pattern can not be used to prune the
      // expanded pattern
      continue;
    }
    assert(!GUNDAM::IsPath<true>(expand_node.pattern()));
    // there is no-way for this pattern to find legal rhs literal
    // prune it from expanded patterns
    std::vector<DfsCodeType> dfs_code_set_for_pattern 
             = GUNDAM::GetDFSCode(expand_node.pattern());
    for (const auto& dfs_code : dfs_code_set_for_pattern) {
      // can successfully erase or not 
      dfs_code_map.erase(dfs_code);
    }
  }
  util::Info("dfs code for pruning generated");

  // ###########################################################
  // ##                                                       ##
  // ##  dfs_code_map contains a map for each dfscode to       ##
  // ##  a set of idx of expanded pattern. Since each         ##
  // ##  expanded pattern can generate multiple dfscodes,     ##
  // ##  the same idx can be contained in different dfscodes. ##
  // ##  Use disjoint set to merge the set of dfscode that    ##
  // ##  contains same idx                                    ##
  // ##                                                       ##
  // ##    should merge all nodes with same dfs_code          ##
  // ##    together here from:                                ##
  // ##       dfs code0 : [1, 2, 3]                           ##
  // ##       dfs code1 : [1, 4, 5]                           ##
  // ##       merge to [1, 2, 3, 4, 5]                        ##
  // ##                                                       ##
  // ###########################################################
  GUNDAM::DisjointSet<size_t> disjoint_set(expand_level.size());
  // since some expanded pattern might have been pruned, needs to 
  // mark the set of patterns that are remained
  std::vector<bool> expand_node_considered(expand_level.size(), false);
  util::Info("begin to merge same dfs code");
  for (const auto& [dfscode, expand_node_id_set]
                  : dfs_code_map) {
    for (const auto& expand_node_id
                   : expand_node_id_set) {
      expand_node_considered[expand_node_id] = true;
      disjoint_set.Merge(*expand_node_id_set.begin(),
                          expand_node_id);
    }
  }
  util::Info("same dfs code merged");
  // free-up space
  dfs_code_map.clear();

  // ########################################################
  // ##  based on the generated disjoint_set, collect the  ##
  // ##  set of expanded pattern to merge                  ##
  // ########################################################


  // use vector at the beginning
  std::map<size_t, // representative idx
           std::vector<size_t>> // set of idx to merge
       merged_expand_node_id_set;
  util::Info("begin to generate the set of nodes to merge");
  for (size_t expand_node_id = 0;
              expand_node_id < expand_level.size();
              expand_node_id++) {
    if (!expand_node_considered[expand_node_id]){
      // this expanded node does not need to be considered
      continue;
    }
    // would an empty std::vector<size_t> if dfs_code 
    // is not contained in dfs_code_map
    merged_expand_node_id_set[disjoint_set.Find(expand_node_id)]
                                  .emplace_back(expand_node_id);
  }
  util::Info("set of nodes to merge generated");

  // for openmp
  std::vector<std::pair<size_t,  // representative idx
            std::vector<size_t>>> // set of idx to merge
       merged_expand_node_id_vector_set;
  merged_expand_node_id_vector_set.reserve(merged_expand_node_id_set.size());
  for (auto& [ representative_id, expand_node_id_group ]
                         : merged_expand_node_id_set) {
    merged_expand_node_id_vector_set.emplace_back(std::move(representative_id),
                                                  std::move(expand_node_id_group));
  }

  // ##############################################
  // ##  generate the next level merged pattern  ##
  // ##############################################
  ExpandTreeLevelType temp_expand_level;
  temp_expand_level.Reserve(merged_expand_node_id_vector_set.size());
  util::Info("begin to par merge");

  omp_lock_t     add_expand_tree_node_lock; 
  omp_init_lock(&add_expand_tree_node_lock);

  #pragma omp parallel for schedule(static) 
  for (size_t merged_expand_node_id_vector_set_idx = 0;
              merged_expand_node_id_vector_set_idx
            < merged_expand_node_id_vector_set.size();
              merged_expand_node_id_vector_set_idx++) {
    const auto& representative_id       = merged_expand_node_id_vector_set[merged_expand_node_id_vector_set_idx].first;
    const auto&    expand_node_id_group = merged_expand_node_id_vector_set[merged_expand_node_id_vector_set_idx].second;

    const GraphPatternType& representative_graph_patten_ref
        = expand_level.node(representative_id).pattern();
        
    #ifndef NDEBUG
    util::Debug("expand_node_id_group.size(): "
               + std::to_string(expand_node_id_group.size()));
    util::Debug("disjoint_set.Find(expand_node_id): "
               + std::to_string(disjoint_set.Find(representative_id)));

    assert(!expand_node_id_group.empty());

    // expand_node in the same group should have the same pattern
    for (const auto& expand_node_id : expand_node_id_group) {
      std::string str;
      util::Debug("pattern0: ");
      str << expand_level.node(expand_node_id).pattern();
      util::Debug(str);

      util::Debug("pattern1: ");
      str.clear();
      str << representative_graph_patten_ref;
      util::Debug(str);

      util::Debug("dfs_code_set0: ");
      const auto dfs_code_set0 = GUNDAM::GetDFSCode(expand_level.node(expand_node_id).pattern());
      for (const auto& dfs_code : dfs_code_set0) {
        util::Debug("\t" + dfs_code.ToString());
      }
      util::Debug("dfs_code_set1: ");
      const auto dfs_code_set1 = GUNDAM::GetDFSCode(representative_graph_patten_ref);
      for (const auto& dfs_code : dfs_code_set1) {
        util::Debug("\t" + dfs_code.ToString());
      }
      assert(GUNDAM::SamePattern(expand_level.node(expand_node_id).pattern(),
                                 representative_graph_patten_ref));
    }
    #endif // NDEBUG

    { // to verify whether the links are legally generated from two side
      std::vector<GraphPatternType> 
             connected_components(GUNDAM::ConnectedComponent<true>(representative_graph_patten_ref));
      std::vector<size_t> can_expand_from(connected_components.size(), 
                                                  expand_level.size() + 1); // intialized as "too many"

      std::vector<GraphPatternType> remove_vertex_graph_pattern_set;
      for (size_t cc_idx = 0; 
                  cc_idx < connected_components.size(); 
                  cc_idx++) {
        const auto& cc_ref = connected_components[cc_idx];
        if (cc_ref.CountVertex() == 1) {
          // has only one vertex
          continue;
        }
        // has more than one vertex, to find whether this is 
        // path or star
        const auto [end_vertex_handle_set,
                central_vertex_handle] = GUNDAM::StarEndPoints<true>(cc_ref);
        if (end_vertex_handle_set.empty()) {
          // is not star or path, cannot be guaranteed 
          // that dfs-code can remove all duplicated patterns
          continue;
        }
        assert(end_vertex_handle_set.size() > 0);
        remove_vertex_graph_pattern_set.reserve(remove_vertex_graph_pattern_set.size()
                                                        + end_vertex_handle_set.size());
        for (const auto& end_vertex_handle 
                       : end_vertex_handle_set) {
          remove_vertex_graph_pattern_set.emplace_back(representative_graph_patten_ref);
          remove_vertex_graph_pattern_set.back()
                                         .FindVertex( end_vertex_handle->id());
          remove_vertex_graph_pattern_set.back()
                                         .EraseVertex(end_vertex_handle->id());
        }
        assert(remove_vertex_graph_pattern_set.size() >= 1);
        assert(remove_vertex_graph_pattern_set.size() <= end_vertex_handle_set.size());
      }

      if (!remove_vertex_graph_pattern_set.empty()) {
        GUNDAM::DeduplicatePatterns(remove_vertex_graph_pattern_set);

        std::set<size_t> expand_vertex_id_set;
        for (const auto& expand_node_id : expand_node_id_group) {
          expand_vertex_id_set.emplace(expand_from_pattern[expand_node_id]);
        }

        // omp_set_lock(&add_expand_tree_node_lock);
        // std::string representative_graph_patten_ref_str;
        // representative_graph_patten_ref_str << representative_graph_patten_ref;
        // util::Info("representative_graph_patten_ref_str: " + representative_graph_patten_ref_str);
        // util::Info("expand_vertex_id_set.size(): " + std::to_string(expand_vertex_id_set.size()));
        // util::Info("can_expand_from_count:       " + std::to_string(remove_vertex_graph_pattern_set.size()));
        // if (expand_vertex_id_set.size() > remove_vertex_graph_pattern_set.size()) {
        //   std::vector<DfsCodeType> representative_graph_patten_ref_dfs_code_set
        //       = GUNDAM::GetDFSCode(representative_graph_patten_ref);
        //   util::Info("representative_graph_patten_ref_dfs_code_set: ");
        //   for (const auto& dfs_code : representative_graph_patten_ref_dfs_code_set) {
        //     util::Info("\t" + dfs_code.ToString());
        //   }
        //   for (const auto& expand_node_id : expand_node_id_group) {
        //     std::string expand_node_pattern_str;
        //     expand_node_pattern_str << expand_level.node(expand_node_id).pattern();
        //     util::Info("\texpand_node_pattern_str: " + expand_node_pattern_str);

        //     std::vector<DfsCodeType> expand_level_node_pattern_dfs_code_set
        //         = GUNDAM::GetDFSCode(expand_level.node(expand_node_id).pattern());

        //     util::Info("\texpand_level_node_pattern_dfs_code_set: ");
        //     for (const auto& dfs_code : expand_level_node_pattern_dfs_code_set) {
        //       util::Info("\t\t" + dfs_code.ToString());
        //     }

        //     std::string pattern_expand_from_str;
        //     pattern_expand_from_str << before_expand_level.node(
        //                                       expand_from_pattern[expand_node_id]).pattern();

        //     util::Info("\tpattern_expand_from id: " 
        //          + std::to_string(expand_from_pattern[expand_node_id]));

        //     util::Info("\tpattern_expand_from_str: " 
        //                 + pattern_expand_from_str);

        //     std::vector<DfsCodeType> pattern_expand_pattern_dfs_code_set
        //         = GUNDAM::GetDFSCode( before_expand_level.node(
        //                                      expand_from_pattern[expand_node_id]).pattern());

        //     util::Info("\tpattern_expand_from_dfs_code_set: ");
        //     for (const auto& dfs_code : pattern_expand_pattern_dfs_code_set) {
        //       util::Info("\t\t" + dfs_code.ToString());
        //     }
        //   }
        //   getchar();
        // }
        // omp_unset_lock(&add_expand_tree_node_lock);

        if (expand_vertex_id_set.size() < remove_vertex_graph_pattern_set.size()) {
          continue;
        }
      }
    }

    omp_set_lock(&add_expand_tree_node_lock);
    auto& added_node = temp_expand_level.AddExpandTreeNode(expand_level.node(representative_id));
    omp_unset_lock(&add_expand_tree_node_lock);
    added_node.rhs_literals().clear();

    std::set<LiteralInfoType> added_rhs_literal_set;

    // ##############################################
    // ##  merge the literal from all other nodes  ##
    // ##  onto the expanded pattern               ##
    // ##############################################
    for (const auto& expand_node_id
                   : expand_node_id_group) {
      auto& node_to_merge = expand_level.node(expand_node_id);
      assert(GUNDAM::SamePattern(added_node.pattern(),
                              node_to_merge.pattern()));

      if (node_to_merge.rhs_literals().empty()) {
        // nothing to merge
        continue;
      }

      auto& pattern_before_expand = before_expand_level.node(
                       expand_from_pattern[expand_node_id]).pattern();

      using MatchType = GUNDAM::Match<GraphPatternType, 
                                      GraphPatternType>;

      std::function<bool(const MatchType&)> prune_nothing_callback
        = [](const MatchType& match_state) -> bool {
          // prune nothing, continue the matching
          return false;
        };

      std::function<bool(const MatchType&)> merge_literal_callback
       = [&added_rhs_literal_set,
          &node_to_merge](const MatchType& match_state) -> bool {
          // only consider rhs literal now
          util::Debug("node_to_merge.rhs_literals().size(): "
                      + std::to_string(node_to_merge.rhs_literals().size()));
          for (const auto& rhs_literal : node_to_merge.rhs_literals()) {
            #ifndef NDEBUG
            std::string rhs_literal_str;
            rhs_literal_str << rhs_literal;
            util::Debug("\t*" + rhs_literal_str);
            #endif // NDEBUG
            auto rhs_literal_on_dst
               = rhs_literal.template MapTo<GraphPatternType>(match_state);
            assert(rhs_literal_on_dst.literal_type() != gar::LiteralType::kNoneLiteral);
            #ifndef NDEBUG
            std::string rhs_literal_on_dst_str;
            rhs_literal_on_dst_str << rhs_literal_on_dst;
            util::Debug("\t#" + rhs_literal_on_dst_str);
            #endif // NDEBUG
            added_rhs_literal_set.emplace(rhs_literal_on_dst);
          }
          // continue matching
          return true;
        };

      util::Debug("merge expand node: " + std::to_string(expand_node_id) );
      auto ret = GUNDAM::MatchUsingMatch<
                 GUNDAM::MatchSemantics::kIsomorphism,
                 GUNDAM::MatchAlgorithm::kVf2,
                 GUNDAM::MergeNecConfig::kNotMerge>(
                                  pattern_before_expand,
                                   added_node.pattern(),
                                 prune_nothing_callback,
                                 merge_literal_callback);
      util::Debug("merged added_rhs_literal_set.size(): " 
                 + std::to_string(added_rhs_literal_set.size()));
      assert(ret >= 1);
    }

    if (restriction.gcr()) {
      // ##########################################
      // ##  file out all literals that are not  ##
      // ##  contained in the end vertexes       ##
      // ##########################################
      std::vector<GraphPatternType> connected_components 
        = GUNDAM::ConnectedComponent(added_node.pattern());
      assert(connected_components.size() == 2);
      std::vector<
      std::vector<typename GUNDAM::VertexID<GraphPatternType>::type>> set_of_end_vertexes_id_set;
      for (const auto& cc : connected_components) {
        const auto [end_vertex_handle_set,
                central_vertex_handle] = GUNDAM::StarEndPoints<true>(cc);
        assert(end_vertex_handle_set.size() >= 2);
        assert(end_vertex_handle_set.size() == 2 || central_vertex_handle);
        set_of_end_vertexes_id_set.emplace_back();
        auto&  end_vertexes_id_set = set_of_end_vertexes_id_set.back();
        end_vertexes_id_set.reserve(end_vertex_handle_set.size());
        for (const auto& vertex_handle 
                   : end_vertex_handle_set) {
          end_vertexes_id_set.emplace_back(vertex_handle->id());
        }
        std::sort(end_vertexes_id_set.begin(),
                  end_vertexes_id_set.end());
        assert(std::is_sorted(end_vertexes_id_set.begin(),
                              end_vertexes_id_set.end()));
      }
      assert(set_of_end_vertexes_id_set.size() == 2);
      util::Debug("gcr: before added_rhs_literal_set.size(): "
              + std::to_string(added_rhs_literal_set.size()));
      for (auto added_rhs_literal_it  = added_rhs_literal_set.begin();
                added_rhs_literal_it != added_rhs_literal_set.end();) {

        const auto& added_rhs_literal = *added_rhs_literal_it;

        const auto vertex_id_set = added_rhs_literal.vertex_id_set();

        if (vertex_id_set.size() == 1) {
          // is constant rule, would not require it be added across
          // different stars, preserve this literal
          added_rhs_literal_it++;
          continue;
        }

        assert(vertex_id_set.size() == 2);
        const auto& vertex_id_0 = vertex_id_set[0];
        const auto& vertex_id_1 = vertex_id_set[1];

        bool contained_in_end_points_of_different_cc = false;
        for (size_t cc_idx = 0; cc_idx < 2 /* at most have two star */; 
                    cc_idx++) {
          const auto& end_vertexes_id_set 
             = set_of_end_vertexes_id_set[cc_idx];
          if (std::binary_search(end_vertexes_id_set.begin(),
                                 end_vertexes_id_set.end(),
                                     vertex_id_0)) {
            // to find whether the other literal is also contained
            // in this component
            const auto& the_other_set_of_end_vertexes_id_set 
                                = set_of_end_vertexes_id_set[1 - cc_idx];
            if (std::binary_search(the_other_set_of_end_vertexes_id_set.begin(),
                                   the_other_set_of_end_vertexes_id_set.end(),
                                                        vertex_id_1)) {
              contained_in_end_points_of_different_cc = true;
              break;
            }
            break;
          }
        }
        if (contained_in_end_points_of_different_cc) {
          // preserve this literal
          added_rhs_literal_it++;
          continue;
        }
        added_rhs_literal_it = added_rhs_literal_set.erase(added_rhs_literal_it);
      }
      util::Debug("gcr: after added_rhs_literal_set.size(): "
             + std::to_string(added_rhs_literal_set.size()));
    }

    // ############################################################
    // ##  map the literal from the expanded back to the         ##
    // ##  patterns before expand to find whether the merged     ## 
    // ##  literal exists on ever pattern, can be pruned if not  ## 
    // ##  Can prune the literals exists on one pattern but      ##
    // ##  have already been pruned by other                     ##
    // ############################################################
    if (!added_rhs_literal_set.empty()) {
      for (const auto& expand_node_id
                     : expand_node_id_group) {
        auto& node_to_merge = expand_level.node(expand_node_id);
        assert(GUNDAM::SamePattern(added_node.pattern(),
                                node_to_merge.pattern()));

        auto& pattern_before_expand = before_expand_level.node(
                                             expand_from_pattern[expand_node_id]).pattern();

        using MatchType = GUNDAM::Match<GraphPatternType, 
                                        GraphPatternType>;

        std::function<bool(const MatchType&)> prune_nothing_callback
          = [](const MatchType& match_state) {
            // prune nothing, continue the matching
            return false;
          };

        std::function<bool(const MatchType&)> erase_literal_callback
         = [&added_rhs_literal_set,
            &node_to_merge](const MatchType& match_state) {
            MatchType match_state_reverse = match_state.Reverse();
            // match all added rhs literal back to node_to_merge.pattern()
            for (auto added_rhs_literal_it  = added_rhs_literal_set.begin();
                      added_rhs_literal_it != added_rhs_literal_set. end ();) {
              auto literal_on_node_to_merge
                 = added_rhs_literal_it->template MapTo<GraphPatternType>(match_state_reverse);
              if (literal_on_node_to_merge.literal_type()
                 == gar::LiteralType::kNoneLiteral) {
                // this literal is not matched by the match_state_reverse
                added_rhs_literal_it++;
                continue;
              }
              // this literal is matched by the match_state_reverse
              bool has_same_literal = false;
              for (const auto& node_to_merge_rhs
                             : node_to_merge.rhs_literals()) {
                if (literal_on_node_to_merge
                            == node_to_merge_rhs){
                  has_same_literal = true;
                  break;
                }
              }
              #ifndef NDEBUG
              std::string literal_on_node_to_merge_str;
              literal_on_node_to_merge_str << literal_on_node_to_merge;
              util::Debug("literal_on_node_to_merge_str: "
                         + literal_on_node_to_merge_str);
              util::Debug("has_same_literal: "
                         + std::to_string(has_same_literal));
              for (const auto& node_to_merge_rhs
                             : node_to_merge.rhs_literals()) {
                std::string node_to_merge_rhs_str;
                node_to_merge_rhs_str << node_to_merge_rhs;
                util::Debug("\tnode_to_merge_rhs: "
                             + node_to_merge_rhs_str);
              }
              #endif // NDEBUG
              if (has_same_literal) {
                // this literal can be considered
                added_rhs_literal_it++;
                continue;
              }
              // does not have the same literal,
              // which means it is literal can be pruned by node_to_merge
              added_rhs_literal_it = added_rhs_literal_set.erase(added_rhs_literal_it);
            }
            // continue matching
            return true;
          };

        std::string pattern_before_expand_str;
        pattern_before_expand_str << pattern_before_expand;

        std::string added_node_pattern_str;
        added_node_pattern_str << added_node.pattern();

        util::Debug("pattern_before_expand_str: "
                   + pattern_before_expand_str);

        util::Debug("added_node_pattern_str: "
                   + added_node_pattern_str);

        auto ret = GUNDAM::MatchUsingMatch<
                   GUNDAM::MatchSemantics::kIsomorphism,
                   GUNDAM::MatchAlgorithm::kVf2,
                   GUNDAM::MergeNecConfig::kNotMerge>(
                                    pattern_before_expand,
                                      added_node.pattern(),
                                    prune_nothing_callback,
                                    erase_literal_callback);
      }
      util::Debug("after pruned all non-legal ones added_rhs_literal_set.size(): "
                                  + std::to_string(added_rhs_literal_set.size()));
      // has merged all literals and pruned all non-legal ones
      // add those literals back to the added expand node
      for (const auto& rhs_literal_info
               : added_rhs_literal_set) {
        added_node.AddRhsLiteral(rhs_literal_info);
      }
    }

    if (restriction.horn_rule()) {

      added_node.ClearRhsLiteral();
    }

    util::Debug("added_node.const_rhs_literals().size(): "
                + std::to_string(added_node.const_rhs_literals().size()));
    #ifndef NDEBUG
    for (const auto& rhs_literal : added_node.const_rhs_literals()) {
      std::string rhs_literal_str;
      rhs_literal_str << rhs_literal;
      util::Debug("\t" + rhs_literal_str);
    }
    #endif // NDEBUG
    // #########################################################
    // ##  add new literal on the end points if it is a link  ##
    // #########################################################
    auto [ link_src_handle,
           link_dst_handle ] = GUNDAM::PathEndPoints<true>(added_node.pattern());
    if (!link_src_handle) {
      // the expanded pattern is not link
      assert(!link_dst_handle);
      // no literal need to be added
      continue;
    }
    // util::Info("is link");
    // util::Info("<src_id, dst_id>: <" + std::to_string(link_src_handle->id())
    //                            + "," + std::to_string(link_dst_handle->id())
    //                            + ">");
    // add literal bewteen the two end points
    std::vector<GraphPatternVertexHandle> vertex_pair;
    vertex_pair.reserve(2);
    vertex_pair.emplace_back(link_src_handle);
    vertex_pair.emplace_back(link_dst_handle);
    std::vector<LiteralInfoType> literals_between_two_end_points;
    ExpandRhsLiteral(added_node,
                     vertex_pair,
                     graph_basic_statistics,
                     ml_edge_label_set, 
                     ml_edge_type_set,
                     ml_to_normal_edge_label,
                        restriction,
                    literals_between_two_end_points);
    util::Debug("literals_between_two_end_points.size(): "
                + std::to_string(literals_between_two_end_points.size()));
    for (auto literal_it  = literals_between_two_end_points.begin();
              literal_it != literals_between_two_end_points.end();) {
      if (literal_it->literal_type() != gar::LiteralType::kEdgeLiteral) {
        // only edge literal can duplicate when a new edge that has been
        // added into the pattern
        assert(added_node.pattern().CountVertex() == 2
            && added_node.pattern().CountEdge()   == 1);
        literal_it++;
        continue;
      }
      if (GUNDAM::HasEdge<const GraphPatternType>(
                    added_node.const_pattern().FindVertex(literal_it->x_id()),
                    literal_it->edge_label(),
                    added_node.const_pattern().FindVertex(literal_it->y_id()))) {
        literal_it = literals_between_two_end_points.erase(literal_it);
        continue;
      }
      literal_it++;
    }
    util::Debug("literals_between_two_end_points.size() after remove duplicate: "
                + std::to_string(literals_between_two_end_points.size()));
    for (const auto& literal_info 
                   : literals_between_two_end_points) {
      added_node.AddRhsLiteral(literal_info);
    }
    util::Debug("literals_between_two_end_points: ");
    for (const auto& rhs_literal : literals_between_two_end_points) {
      std::string rhs_literal_str;
      rhs_literal_str << rhs_literal;
      util::Debug("\t" + rhs_literal_str);
    }
    std::string pattern_str;
    pattern_str << added_node.pattern();
    util::Debug("pattern_str: " + pattern_str);
    for (const auto& rhs_literal : added_node.const_rhs_literals()) {
      std::string rhs_literal_str;
      rhs_literal_str << rhs_literal;
      util::Debug("\t" + rhs_literal_str);
    }
  }
  util::Info("merged");
  std::swap(expand_level, temp_expand_level);
  return;
}

} // namespace _gar_discover

} // namespace grape

#endif // _UNIQUE_PATTERN_H
