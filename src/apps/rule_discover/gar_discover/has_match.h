#ifndef EXAMPLES_ANALYTICAL_HAS_MATCH_
#define EXAMPLES_ANALYTICAL_HAS_MATCH_

#include <algorithm>

#include "gar/gar_supp.h"

#include "gundam/tool/connected_component.h"
#include "gundam/tool/sub_graph_of.h"
#include "gundam/tool/has_edge.h"

namespace grape {

namespace _gar_discover {

namespace _has_match {
static constexpr float does_not_return_confidence = 0.0;
};

template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
std::pair<uint64_t, float> HasMatch(GraphPatternType& pattern,
 const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literal_set,
             const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal,
                                       DataGraphType& data_graph,
       const std::map<typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type,
                      typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type>& normal_to_ml_edge_label,
                  const CandidateSetContainer& pattern_candidate_set,
                  const CandidateSetContainer& pattern_candidate_set_removed_nec,
                  const uint64_t supp_bound,
                  const bool only_verify_is_legal,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_statistic,
                  const Restriction<GraphPatternType,
                                       DataGraphType>& restrictions,
                  const double time_limit,
                  const double time_limit_per_supp,
                  const bool store_match,
    std::vector<std::pair<std::vector<typename GUNDAM::VertexID<DataGraphType>::type>,
                          std::vector<typename GUNDAM::EdgeID<DataGraphType>::type>>>& match_id_vec,
                  double& probability) {
                  
  const bool kConsiderConfidenceBound 
          = !only_verify_is_legal
          && restrictions.specified_confidence_bound();

  assert(literal_set.size() >= 1);
  #ifndef NDEBUG 
  bool _has_same = false;
  for (const auto& literal : literal_set) {
    if (literal == rhs_literal) {
      assert(!_has_same); // should only have one same
      _has_same = true;
    }
  }
  assert(_has_same);
  #endif // NDEBUG 
  
  if (restrictions.gcr()) {
    /* ######################################################## *
     * ##  pruning:                                          ## *
     * ##    since gcr utilize homomorphism and consider     ## *
     * ##    star-pattern only, when a path without literal  ## *
     * ##    can be contained in another path, this gcr      ## *
     * ##    has been evalauted before and satisfy the       ## *
     * ##    support bound. Need to be evaluated again       ## *
     * ######################################################## */
    // first collect id set of vertexes that are contained in
    // the literals
    std::vector<typename GUNDAM::VertexID<GraphPatternType>::type > 
      vertex_id_set_in_all_literals
                     = rhs_literal.vertex_id_set();
    for (const auto& literal : literal_set) {
      auto vertex_id_set = literal.vertex_id_set();
      vertex_id_set_in_all_literals.insert(
      vertex_id_set_in_all_literals.end(),
        std::make_move_iterator(vertex_id_set.begin()),
        std::make_move_iterator(vertex_id_set.end()));
    }
    std::sort(vertex_id_set_in_all_literals.begin(),
              vertex_id_set_in_all_literals.end());

    vertex_id_set_in_all_literals.erase(
        std::unique(vertex_id_set_in_all_literals.begin(), 
                    vertex_id_set_in_all_literals.end()), 
                    vertex_id_set_in_all_literals.end());

    std::vector<GraphPatternType> 
          connected_components(GUNDAM::ConnectedComponent(pattern));

    for (auto& connected_component : connected_components) {
      const auto [end_vertex_handle_set,
              central_vertex_handle] = GUNDAM::StarEndPoints<true>(connected_component);

      GUNDAM::Match<GraphPatternType,
                    GraphPatternType> partial_match;

      if (kConsiderConfidenceBound) {
        if (central_vertex_handle) {
          partial_match.AddMap(connected_component.FindVertex(central_vertex_handle->id()),
                               connected_component.FindVertex(central_vertex_handle->id()));
        }
        else {
          assert(end_vertex_handle_set.size() == 2);
        }
      }

      for (const auto& vertex_id : vertex_id_set_in_all_literals) {
        if (!connected_component.FindVertex(vertex_id)) {
          continue;
        }
        partial_match.AddMap(connected_component.FindVertex(vertex_id),
                             connected_component.FindVertex(vertex_id));
      }

      auto ret = GUNDAM::MatchUsingMatch<
                 GUNDAM::MatchSemantics::kHomomorphism>(connected_component,
                                                        connected_component,
                                                        partial_match);

      assert(ret >= 1);

      if (ret >= 2) {
        // has been considered before, satisfy support bound but not
        // confidence bound (if specified)
        return std::pair(supp_bound, _has_match::does_not_return_confidence);
      }
    }
  }

  gar::GraphAssociationRule<GraphPatternType,
                                DataGraphType> temp_gar(pattern);
  temp_gar.AddY(rhs_literal);

  const auto kRhsLiteralStandAloneInfo
          = (*temp_gar.y_literal_set()
                      .begin())->stand_alone_info();

  CandidateSetContainer temp_candidate_set,
                        temp_candidate_set_removed_nec;

  auto temp_candidate_set_ptr 
     = std::addressof(pattern_candidate_set),
       temp_candidate_set_removed_nec_ptr 
     = std::addressof(pattern_candidate_set_removed_nec);

  if (!restrictions.specified_rhs_literal_set().empty()) {
    // specified rhs literal set
    if (restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
      assert(kRhsLiteralStandAloneInfo.VertexNum() == 2);
      const auto& candidate_set_for_vertex_in_rhs_literal 
                = restrictions
                 .candidate_set_for_vertex_in_rhs_literal(kRhsLiteralStandAloneInfo);
                
      temp_candidate_set             = pattern_candidate_set;
      temp_candidate_set_removed_nec = pattern_candidate_set_removed_nec;
      
      const auto x_handle = pattern.FindVertex(rhs_literal.x_id());
      const auto y_handle = pattern.FindVertex(rhs_literal.y_id());

      std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type> 
             literal_vertex_handle_set{x_handle, y_handle};

      assert(literal_vertex_handle_set.size() 
          == candidate_set_for_vertex_in_rhs_literal.size());

      // for each vertex contained in rhs literal, intersect its
      // candidate set with the vertex set apeared in the specified
      // support set
      for (size_t vertex_idx = 0; 
                  vertex_idx < literal_vertex_handle_set.size(); 
                  vertex_idx++) {
        const auto& vertex_handle = literal_vertex_handle_set[vertex_idx];
        const auto& specified_candidate_set_for_vertex 
                            = candidate_set_for_vertex_in_rhs_literal[vertex_idx];

         auto  vertex_handle_candidate_set_it  = temp_candidate_set.find(vertex_handle);
        assert(vertex_handle_candidate_set_it != temp_candidate_set.end());
        auto& vertex_handle_candidate_set 
            = vertex_handle_candidate_set_it->second;

         auto  vertex_handle_candidate_set_removed_nec_it  = temp_candidate_set_removed_nec.find(vertex_handle);
        assert(vertex_handle_candidate_set_removed_nec_it != temp_candidate_set_removed_nec.end());
        auto& vertex_handle_candidate_set_removed_nec 
            = vertex_handle_candidate_set_removed_nec_it->second;

        std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> intersetct_candidate;
        intersetct_candidate.reserve(
                  std::min(vertex_handle_candidate_set.size(),
                    specified_candidate_set_for_vertex.size()));
        std::set_intersection(vertex_handle_candidate_set.begin(),
                              vertex_handle_candidate_set. end (),
                       specified_candidate_set_for_vertex.begin(),
                       specified_candidate_set_for_vertex. end (),
                              std::back_inserter(intersetct_candidate));
        intersetct_candidate.swap(vertex_handle_candidate_set);
        intersetct_candidate.clear();

        intersetct_candidate.reserve(
                  std::min(vertex_handle_candidate_set_removed_nec.size(),
                               specified_candidate_set_for_vertex.size()));
        std::set_intersection(vertex_handle_candidate_set_removed_nec.begin(),
                              vertex_handle_candidate_set_removed_nec. end (),
                                  specified_candidate_set_for_vertex.begin(),
                                  specified_candidate_set_for_vertex. end (),
                              std::back_inserter(intersetct_candidate));
        intersetct_candidate.swap(vertex_handle_candidate_set_removed_nec);
        intersetct_candidate.clear();
      }
      temp_candidate_set_ptr             = std::addressof(temp_candidate_set);
      temp_candidate_set_removed_nec_ptr = std::addressof(temp_candidate_set_removed_nec);
    }
  }

  assert(temp_candidate_set_ptr);
  assert(temp_candidate_set_removed_nec_ptr);

  const CandidateSetContainer& processed_pattern_candidate_set 
                                         = *temp_candidate_set_ptr,
                               processed_pattern_candidate_set_removed_nec
                                         = *temp_candidate_set_removed_nec_ptr;

  if (restrictions.gcr()) {
    // when has only one non-constant literal, can evaluate
    // it through the simulation
    size_t lhs_non_constant_literal_counter = 0,
               lhs_constant_literal_counter = 0,
               lhs_variable_literal_counter = 0;
    gar::LiteralInfo<GraphPatternType,
                        DataGraphType> non_constant_literal;
    assert(non_constant_literal.literal_type() == gar::LiteralType::kNoneLiteral);
    for (const auto& literal : literal_set) {
      assert(literal.literal_type() != gar::LiteralType::kNoneLiteral);
      if (literal.literal_type() != gar::LiteralType::kConstantLiteral) {
        non_constant_literal = literal;
      }
      /* ################################## *
       * ##  statistic the lhs literals  ## *
       * ################################## */
      if (kConsiderConfidenceBound) {
        if (literal == rhs_literal) {
          // does not add rhs literal to lhs
          continue;
        }
      }
      if (literal.literal_type() == gar::LiteralType::kConstantLiteral) {
        lhs_constant_literal_counter++;
        continue;
      }
      lhs_non_constant_literal_counter++;
      if (literal.literal_type() == gar::LiteralType::kVariableLiteral) {
        lhs_variable_literal_counter++;
        continue;
      }
    }

    if ( lhs_non_constant_literal_counter == 0
     || (lhs_non_constant_literal_counter == 1
     && !kConsiderConfidenceBound)
  /* || (lhs_non_constant_literal_counter == 1 // only support one variable literal in lhs now 
          && lhs_variable_literal_counter == 1) */) {
      // can utilize simulation
      const CandidateSetContainer* filted_constant_candidate_set_ptr
       = std::addressof(processed_pattern_candidate_set);

      CandidateSetContainer candidate_set_satisfy_constant_rules;
      if (lhs_constant_literal_counter > 0) {
        // has constant_literal in lhs, filt 
        candidate_set_satisfy_constant_rules = processed_pattern_candidate_set;
        for (const auto& literal : literal_set) {
          if (literal.literal_type() != gar::LiteralType::kConstantLiteral) {
            continue;
          }
          if (kConsiderConfidenceBound) {
            // does not add rhs literal to lhs
            if (literal == rhs_literal) {
              continue;
            }
          }
          // filt out all candidate vertexes not satisfy the constant rules
          const auto x_handle = pattern.FindVertex(literal.x_id());
          const auto x_attr_key   = literal.x_attr_key();
          const auto x_value_str  = literal.c_str();
          const auto x_value_type = literal.data_type();
           auto  x_candidate_set_it  = candidate_set_satisfy_constant_rules.find(x_handle);
          assert(x_candidate_set_it != candidate_set_satisfy_constant_rules.end());
          auto& x_candidate_set 
              = x_candidate_set_it->second;
          std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> satisfy_literal_x_candidate_set;
          satisfy_literal_x_candidate_set.reserve(
                          x_candidate_set.size());
          for (const auto& dst_handle : x_candidate_set) {
            const auto x_dst_attr_handle = dst_handle->FindAttribute(x_attr_key);
            if (!x_dst_attr_handle) {
              // does not has attribute, not satisfy this literal
              continue;
            }
            // has this attribute
            if (x_dst_attr_handle->value_type() != x_value_type
             || x_dst_attr_handle->value_str()  != x_value_str) {
              // has different type or value, not satisfy the literal
              continue;
            }
            // satisfy literal
            satisfy_literal_x_candidate_set.emplace_back(dst_handle);
          }
          x_candidate_set.swap(satisfy_literal_x_candidate_set);
        }
        GUNDAM::_dp_iso_using_match
              ::RefineCandidateSet(pattern, data_graph, 
                                   candidate_set_satisfy_constant_rules);
        filted_constant_candidate_set_ptr = std::addressof(candidate_set_satisfy_constant_rules);
      }
      const CandidateSetContainer& filted_constant_candidate_set
                                = *filted_constant_candidate_set_ptr;
      
      if (rhs_literal.literal_type() == gar::LiteralType::kConstantLiteral) {
        // all literals are constant literal
        assert(non_constant_literal.literal_type() == gar::LiteralType::kNoneLiteral);
        assert(         rhs_literal.literal_type() == gar::LiteralType::kConstantLiteral);
        const auto x_handle = pattern.FindVertex(rhs_literal.x_id());
        assert(x_handle);
         auto  x_candidate_set_it  = filted_constant_candidate_set.find(x_handle);
        assert(x_candidate_set_it != filted_constant_candidate_set.end());
         auto& x_candidate_set 
             = x_candidate_set_it->second;
        if (!kConsiderConfidenceBound) {
          return std::pair(std::min(x_candidate_set.size(), supp_bound),
                    _has_match::does_not_return_confidence);
        }
        if (x_candidate_set.empty()) {
          return std::pair(0, 0);
        }
         auto  x_supp_counter = x_candidate_set.size();
        assert(x_supp_counter > 0);
        size_t xy_supp_counter = 0;
        const auto x_attr_key   = rhs_literal.x_attr_key();
        const auto x_value_str  = rhs_literal.c_str();
        const auto x_value_type = rhs_literal.data_type();
        for (const auto& dst_handle : x_candidate_set) {
          const auto x_dst_attr_handle = dst_handle->FindAttribute(x_attr_key);
          if (!x_dst_attr_handle) {
            // does not has attribute, not satisfy this literal
            continue;
          }
          if (x_dst_attr_handle->value_type() != x_value_type
           || x_dst_attr_handle->value_str()  != x_value_str) {
            // has different type or value, not satisfy the rhs_literal
            continue;
          }
          // satisfy the rhs_literal
          xy_supp_counter++;
        }
        assert(x_supp_counter > 0);
        return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
                                        / ((float)  x_supp_counter));
      }
      assert(non_constant_literal == rhs_literal);
      switch (rhs_literal.literal_type()) {
       case gar::LiteralType::kVariableLiteral: {
         auto  x_handle = pattern.FindVertex(rhs_literal.x_id());
        assert(x_handle);
         auto  y_handle = pattern.FindVertex(rhs_literal.y_id());
        assert(y_handle);

        const auto& x_attr_key = rhs_literal.x_attr_key();
        const auto& y_attr_key = rhs_literal.y_attr_key();

        assert( filted_constant_candidate_set.find(x_handle)
             != filted_constant_candidate_set.end() );
        const auto&  x_candidate_set 
              = filted_constant_candidate_set.find(x_handle)->second;

        assert( filted_constant_candidate_set.find(y_handle)
             != filted_constant_candidate_set.end() );
        const auto&  y_candidate_set 
              = filted_constant_candidate_set.find(y_handle)->second;

        uint64_t x_supp_counter = x_candidate_set.size()
                                * y_candidate_set.size();
        if (x_supp_counter == 0) {
          return std::pair(0, 0);
        }

        std::map<std::string, typename DataGraphType
                                 ::VertexCounterType> x_histogram,
                                                      y_histogram;

        for (const auto& x_candidate : x_candidate_set) {
          const auto x_attr_handle = x_candidate->FindAttribute(x_attr_key);
          if (!x_attr_handle) {
            // does not have this literal
            continue;
          }
          x_histogram[x_attr_handle->value_str()]++;
        }
        for (const auto& y_candidate : y_candidate_set) {
          const auto y_attr_handle = y_candidate->FindAttribute(y_attr_key);
          if (!y_attr_handle) {
            // does not have this literal
            continue;
          }
          y_histogram[y_attr_handle->value_str()]++;
        }

        uint64_t xy_supp_counter = 0;
        for (const auto& [x_str, x_counter] : x_histogram) {
          auto y_histogram_it =  y_histogram.find(x_str);
          if ( y_histogram_it == y_histogram.end() ) {
            continue;
          }
          xy_supp_counter += x_counter * y_histogram_it->second;
          if (kConsiderConfidenceBound) {
            // does not end when reach support bound
            continue;
          }
          if ( xy_supp_counter >= supp_bound ) {
            // has reached the support bound
            assert(!kConsiderConfidenceBound);
            return std::pair(supp_bound, _has_match::does_not_return_confidence);
          }
        }
        if (!kConsiderConfidenceBound) {
          assert(xy_supp_counter < supp_bound);
          return std::pair(xy_supp_counter, _has_match::does_not_return_confidence);
        }
        assert(x_supp_counter > 0);
        return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
                                        / ((float)  x_supp_counter));
       }
       case gar::LiteralType::kEdgeLiteral: {
         auto  x_handle = pattern.FindVertex(rhs_literal.x_id());
        assert(x_handle);
         auto  y_handle = pattern.FindVertex(rhs_literal.y_id());
        assert(y_handle);

        assert( filted_constant_candidate_set.find(x_handle)
             != filted_constant_candidate_set.end() );
        const auto& x_vertex_x_supp_candidate_set 
                           = filted_constant_candidate_set.find(x_handle)->second;

        assert( filted_constant_candidate_set.find(y_handle)
             != filted_constant_candidate_set.end() );
        const auto& y_vertex_x_supp_candidate_set 
                           = filted_constant_candidate_set.find(y_handle)->second;

        uint64_t x_supp_counter = 0;
        if (kConsiderConfidenceBound) {
          // other wise, does not need to process x_supp_counter
          if (!restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
            x_supp_counter = x_vertex_x_supp_candidate_set.size()
                           * y_vertex_x_supp_candidate_set.size();
          }
          else {
            for (const auto& support : restrictions.support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
              assert(support.size() == 2);
              if (!std::binary_search(x_vertex_x_supp_candidate_set.begin(),
                                      x_vertex_x_supp_candidate_set.end(),
                                      support[0])) {
                continue;
              }
              if (!std::binary_search(y_vertex_x_supp_candidate_set.begin(),
                                      y_vertex_x_supp_candidate_set.end(),
                                      support[1])) {
                continue;
              }
              x_supp_counter++;
            }
          }
        }

        uint64_t xy_supp_counter = 0;

        if (!restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
          for (auto vertex_it = data_graph.VertexBegin(x_handle->label());
                   !vertex_it.IsDone();
                    vertex_it++) {
            typename GUNDAM::VertexHandle<DataGraphType>::type
                      vertex_handle = vertex_it;
            if (!std::binary_search(x_vertex_x_supp_candidate_set.begin(),
                                    x_vertex_x_supp_candidate_set.end(),
                                      vertex_handle)) {
              continue;
            }
            for (auto out_edge_it = vertex_it->OutEdgeBegin();
                     !out_edge_it.IsDone();
                      out_edge_it++) {
              if (out_edge_it->label() != rhs_literal.edge_label()) {
                continue;
              }
              if (out_edge_it->dst_handle()->label() != y_handle->label()) {
                continue;
              }
              if (!std::binary_search(y_vertex_x_supp_candidate_set.begin(),
                                      y_vertex_x_supp_candidate_set.end(),
                                      out_edge_it->dst_handle())) {
                continue;
              }
              xy_supp_counter++;
            }
          }
        }
        else {
          for (const auto& support : restrictions.support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
            assert(support.size() == 2);
            if (!std::binary_search(x_vertex_x_supp_candidate_set.begin(),
                                    x_vertex_x_supp_candidate_set.end(),
                                    support[0])) {
              continue;
            }
            if (!std::binary_search(y_vertex_x_supp_candidate_set.begin(),
                                    y_vertex_x_supp_candidate_set.end(),
                                    support[1])) {
              continue;
            }
            if (!GUNDAM::HasEdge<DataGraphType>(
                                support[0],
                                rhs_literal.edge_label(),
                                support[1])) {
              continue;
            }
            xy_supp_counter++;
          }
        }

        if (!kConsiderConfidenceBound) {
          return std::pair(std::min(xy_supp_counter, supp_bound), 
                         _has_match::does_not_return_confidence);
        }
        if (x_supp_counter == 0) {
          return std::pair(xy_supp_counter, 0);
        }
        assert(x_supp_counter > 0);
        return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
                                        / ((float)  x_supp_counter));
       }
       default:
        assert(false);
        break;
      }
      assert(false);
      return std::pair(0, _has_match::does_not_return_confidence);
    }
  }

  // optimization for edge_literal
  bool gar_pattern_added_edge = false;
  GraphPatternType pattern_with_edge_literal(pattern);
  // there is an edge literal in literal_set, then add it as an edge
  // into the pattern for gar
  typename GUNDAM::EdgeID<GraphPatternType>::type 
    new_edge_id = GUNDAM::MaxEdgeId(pattern_with_edge_literal);
  new_edge_id++;
  for (const auto& literal : literal_set) {
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral) {
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               literal.edge_label(),
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
    if (literal.literal_type() == gar::LiteralType::kMlLiteral) {
       auto  ml_edge_label_it  = normal_to_ml_edge_label.find(literal.edge_label());
      assert(ml_edge_label_it != normal_to_ml_edge_label.end());
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               ml_edge_label_it->second,
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
  }
  gar::GraphAssociationRule<GraphPatternType,
                               DataGraphType> gar(pattern_with_edge_literal);
  gar.AddY(rhs_literal);
  // add both x and y literal into x
  for (const auto& literal : literal_set) {
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral){
      // edge literal has already been added into the
      // pattern
      continue;
    }
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    gar.AddX(literal);
  }
  auto& gar_pattern = gar.pattern();
  CandidateSetContainer gar_pattern_candidate_set;
  for (const auto& processed_pattern_candidate 
                 : processed_pattern_candidate_set) {
    auto pattern_vertex_handle = processed_pattern_candidate.first;
    // should not be nullptr
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle 
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    // should not be nullptr
    assert(gar_pattern_vertex_handle);
    auto [ gar_pattern_candidate_set_it,
           gar_pattern_candidate_set_ret ]
         = gar_pattern_candidate_set.emplace(gar_pattern_vertex_handle,
                                       processed_pattern_candidate.second);
    // should added successfully
    assert(gar_pattern_candidate_set_ret);
  }
  assert( processed_pattern_candidate_set.size()
            ==  gar_pattern_candidate_set.size() );

  CandidateSetContainer gar_pattern_candidate_set_removed_nec;
  for (const auto& processed_pattern_candidate 
                 : processed_pattern_candidate_set_removed_nec) {
    auto pattern_vertex_handle = processed_pattern_candidate.first;
    // should not be nullptr
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle 
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    // should not be nullptr
    assert(gar_pattern_vertex_handle);
    auto [ gar_pattern_candidate_set_removed_nec_it,
           gar_pattern_candidate_set_removed_nec_ret ]
         = gar_pattern_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                                   processed_pattern_candidate.second);
    // should added successfully
    assert(gar_pattern_candidate_set_removed_nec_ret);
  }
  assert( processed_pattern_candidate_set_removed_nec.size()
             == gar_pattern_candidate_set_removed_nec.size() );

  if (gar_pattern_added_edge) {
    // further refine the candidate set for more efficient processing
    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                          data_graph, 
                               gar_pattern_candidate_set_removed_nec)) {
      return std::pair(0, _has_match::does_not_return_confidence);
    }
    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                          data_graph, 
                                           gar_pattern_candidate_set)) {
      return std::pair(0, _has_match::does_not_return_confidence);
    }
  }

  uint64_t x_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter,
           &supp_bound,
           &kConsiderConfidenceBound](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || x_supp_counter < supp_bound);
    x_supp_counter++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (x_supp_counter == supp_bound) {
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  uint64_t xy_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter,
           &supp_bound,
           &kConsiderConfidenceBound](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || xy_supp_counter < supp_bound);
    xy_supp_counter++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (xy_supp_counter == supp_bound) {
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  assert(kRhsLiteralStandAloneInfo
     == (*gar.y_literal_set()
             .begin())->stand_alone_info());

  if (!restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
    GUNDAM::Match<GraphPatternType,
                    DataGraphType> match_state;
    if (restrictions.gcr()) {
      // gcr, use homomorphism
      auto [x_supp, xy_supp] = gar::GarSupp<GUNDAM::MatchSemantics::kHomomorphism>(
                                            gar, data_graph, match_state,
                                            gar_pattern_candidate_set,
                                            gar_pattern_candidate_set_removed_nec,
                                            satisfy_x_supp_callback, 
                                            satisfy_xy_supp_callback, 
                                            time_limit,
                                            time_limit_per_supp);
      assert(xy_supp == xy_supp_counter);
      assert( x_supp ==  x_supp_counter);
      if (!kConsiderConfidenceBound) {
        assert( x_supp_counter == xy_supp_counter);
        return std::pair(xy_supp_counter, _has_match::does_not_return_confidence);
      }
      if (x_supp_counter == 0) {
        return std::pair(xy_supp_counter, 0);
      }
      assert(x_supp_counter > 0);
      return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
                                      / ((float)  x_supp_counter));
    }
    // gar or horn rule, use isomorphism
    std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>> match_vec;
    auto [x_supp, xy_supp] = gar::GarSupp<GUNDAM::MatchSemantics::kIsomorphism>(
                                          gar, data_graph, match_state,
                                          gar_pattern_candidate_set,
                                          gar_pattern_candidate_set_removed_nec,
                                          store_match, match_vec, probability,
                                          satisfy_x_supp_callback, 
                                          satisfy_xy_supp_callback, 
                                          time_limit,
                                          time_limit_per_supp);

    assert(xy_supp == xy_supp_counter);
    assert( x_supp ==  x_supp_counter);
    if (!kConsiderConfidenceBound) {
      assert( x_supp_counter == xy_supp_counter);
      return std::pair(xy_supp_counter, _has_match::does_not_return_confidence);
    }
    if (x_supp_counter == 0) {
      return std::pair(xy_supp_counter, 0);
    }
    if (store_match) {
      TransferMatchToID(gar, data_graph, match_vec, match_id_vec);
    }
    assert(x_supp_counter > 0);
    return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
                                    / ((float)  x_supp_counter));
  }

  // specified support set 
  const auto& support_set = restrictions.support_set_for_rhs_literal(kRhsLiteralStandAloneInfo);
  assert(!support_set.empty());

  const auto kRhsLiteralInfo
         = (*gar.y_literal_set()
                .begin())->info();
  std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>> match_vec;
  auto begin_time = std::time(NULL);
  for (const auto& support : support_set) {
    assert(kRhsLiteralStandAloneInfo.VertexNum() == support.size());
    assert(kRhsLiteralStandAloneInfo.VertexNum() == 2);

    if (time_limit > 0 
     && time_limit < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }

    GUNDAM::Match<GraphPatternType,
                     DataGraphType> partial_match_state;

    assert(gar.pattern().FindVertex(kRhsLiteralInfo.x_id()));
    assert(gar.pattern().FindVertex(kRhsLiteralInfo.y_id()));
    partial_match_state.AddMap(gar.pattern().FindVertex(kRhsLiteralInfo.x_id()),
                               support[0]);
    partial_match_state.AddMap(gar.pattern().FindVertex(kRhsLiteralInfo.y_id()),
                               support[1]);
    
    if (restrictions.gcr()) {
      // gcr, use homomorphism
      auto [x_supp, xy_supp] = gar::GarSupp<GUNDAM::MatchSemantics::kHomomorphism>(
                                            gar, data_graph, partial_match_state,
                                            gar_pattern_candidate_set,
                                            gar_pattern_candidate_set_removed_nec,
                                            satisfy_x_supp_callback, 
                                            satisfy_xy_supp_callback, 
                                            time_limit_per_supp,
                                            time_limit_per_supp);
      assert(xy_supp <= x_supp);
      assert(xy_supp <= 1 && x_supp <= 1);
      continue;
    }
    // gar or horn rule, use isomorphism
    auto [x_supp, xy_supp] = gar::GarSupp<GUNDAM::MatchSemantics::kIsomorphism>(
                                          gar, data_graph, partial_match_state,
                                          gar_pattern_candidate_set,
                                          gar_pattern_candidate_set_removed_nec,
                                          store_match, match_vec, probability,
                                          satisfy_x_supp_callback, 
                                          satisfy_xy_supp_callback, 
                                          time_limit_per_supp,
                                          time_limit_per_supp);
    assert(xy_supp <= x_supp);
    assert(xy_supp <= 1 && x_supp <= 1);
  }
  if (!kConsiderConfidenceBound) {
    assert( x_supp_counter == xy_supp_counter);
    return std::pair(xy_supp_counter, _has_match::does_not_return_confidence);
  }
  if (x_supp_counter == 0) {
    return std::pair(xy_supp_counter, 0);
  }

  if (store_match) {
    TransferMatchToID(gar, data_graph, match_vec, match_id_vec);
  }
  assert(x_supp_counter > 0);
  return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
                                  / ((float)  x_supp_counter));
}


template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
std::pair<uint64_t, float> HasMatch(GraphPatternType& pattern,
 const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literal_set,
             const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal,
                                       DataGraphType& data_graph,
       const std::map<typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type,
                      typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type>& normal_to_ml_edge_label,
                  const CandidateSetContainer& pattern_candidate_set,
                  const CandidateSetContainer& pattern_candidate_set_removed_nec,
                  const uint64_t supp_bound,
                  const bool only_verify_is_legal,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_statistic,
                  const Restriction<GraphPatternType,
                                       DataGraphType>& restrictions,
                  const double time_limit,
                  const double time_limit_per_supp) {
  bool store_match = false;
  std::vector<
      std::pair<std::vector<typename GUNDAM::VertexID<DataGraphType>::type>,
                std::vector<typename GUNDAM::EdgeID<DataGraphType>::type>>> match_id_vec;
  double probability = 0.0;
  return HasMatch(pattern, literal_set, rhs_literal, data_graph,
                  normal_to_ml_edge_label, pattern_candidate_set, 
                  pattern_candidate_set_removed_nec,
                  supp_bound, only_verify_is_legal, graph_statistic, restrictions,
                  time_limit, time_limit_per_supp, store_match, match_id_vec, probability);
}


template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
auto HasMatch(const GraphPatternType& pattern,
              const std::vector<gar::LiteralInfo<GraphPatternType,
                                                    DataGraphType>>& literal_set,
                          const gar::LiteralInfo<GraphPatternType,
                                                    DataGraphType>&  rhs_literal,
                      DataGraphType& data_graph,
              const std::map<typename DataGraphType::EdgeType::LabelType,
                             typename DataGraphType::EdgeType::LabelType>& normal_to_ml_edge_label,
        const CandidateSetContainer& pattern_candidate_set,
              const uint64_t supp_bound,
              const bool only_verify_is_legal,
              const double time_limit,
              const double time_limit_per_supp) {
  return HasMatch(pattern,
                  literal_set,
                  rhs_literal,
                  data_graph,
                  normal_to_ml_edge_label,
                  pattern_candidate_set,
                  pattern_candidate_set,
                  supp_bound,
                  only_verify_is_legal,
                  time_limit,
                  time_limit_per_supp);
}



template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
std::pair<std::vector<uint64_t>, std::vector<float>> HasMatchForSet(
 GraphPatternType& pattern,
 const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literal_set,
             const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal,
             const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& new_literals_to_check,
                                       DataGraphType& data_graph,
       const std::map<typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type,
                      typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type>& normal_to_ml_edge_label,
                  const CandidateSetContainer& pattern_candidate_set,
                  const CandidateSetContainer& pattern_candidate_set_removed_nec,
                  const uint64_t supp_bound,
                  const bool only_verify_is_legal,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_statistic,
                  const Restriction<GraphPatternType,
                                       DataGraphType>& restrictions,
                  const double time_limit,
                  const double time_limit_per_supp,
                  bool store_match,
    std::vector<std::vector<std::pair<std::vector<typename GUNDAM::VertexID<DataGraphType>::type>,
                          std::vector<typename GUNDAM::EdgeID<DataGraphType>::type>>>>& match_id_vec_for_all) {
                  
  const bool kConsiderConfidenceBound 
          = !only_verify_is_legal
          && restrictions.specified_confidence_bound();
  assert(new_literals_to_check.size() >= 1);
  assert(literal_set.size() >= 1);
  #ifndef NDEBUG 
  bool _has_same = false;
  for (const auto& literal : literal_set) {
    if (literal == rhs_literal) {
      assert(!_has_same); // should only have one same
      _has_same = true;
    }
  }
  assert(_has_same);
  #endif // NDEBUG 
  
  if (restrictions.gcr()) {
    std::cout << "not suitbale for gcr" << std::endl;
    return std::pair(std::vector<uint64_t>(), std::vector<float>());
  }

  gar::GraphAssociationRule<GraphPatternType,
                                DataGraphType> temp_gar(pattern);
  temp_gar.AddY(rhs_literal);

  const auto kRhsLiteralStandAloneInfo
          = (*temp_gar.y_literal_set()
                      .begin())->stand_alone_info();

  CandidateSetContainer temp_candidate_set,
                        temp_candidate_set_removed_nec;

  auto temp_candidate_set_ptr 
     = std::addressof(pattern_candidate_set),
       temp_candidate_set_removed_nec_ptr 
     = std::addressof(pattern_candidate_set_removed_nec);

  if (!restrictions.specified_rhs_literal_set().empty()) {
    // specified rhs literal set
    if (restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
      assert(kRhsLiteralStandAloneInfo.VertexNum() == 2);
      const auto& candidate_set_for_vertex_in_rhs_literal 
                = restrictions
                 .candidate_set_for_vertex_in_rhs_literal(kRhsLiteralStandAloneInfo);
                
      temp_candidate_set             = pattern_candidate_set;
      temp_candidate_set_removed_nec = pattern_candidate_set_removed_nec;
      
      const auto x_handle = pattern.FindVertex(rhs_literal.x_id());
      const auto y_handle = pattern.FindVertex(rhs_literal.y_id());

      std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type> 
             literal_vertex_handle_set{x_handle, y_handle};

      assert(literal_vertex_handle_set.size() 
          == candidate_set_for_vertex_in_rhs_literal.size());

      // for each vertex contained in rhs literal, intersect its
      // candidate set with the vertex set apeared in the specified
      // support set
      for (size_t vertex_idx = 0; 
                  vertex_idx < literal_vertex_handle_set.size(); 
                  vertex_idx++) {
        const auto& vertex_handle = literal_vertex_handle_set[vertex_idx];
        const auto& specified_candidate_set_for_vertex 
                            = candidate_set_for_vertex_in_rhs_literal[vertex_idx];

         auto  vertex_handle_candidate_set_it  = temp_candidate_set.find(vertex_handle);
        assert(vertex_handle_candidate_set_it != temp_candidate_set.end());
        auto& vertex_handle_candidate_set 
            = vertex_handle_candidate_set_it->second;

         auto  vertex_handle_candidate_set_removed_nec_it  = temp_candidate_set_removed_nec.find(vertex_handle);
        assert(vertex_handle_candidate_set_removed_nec_it != temp_candidate_set_removed_nec.end());
        auto& vertex_handle_candidate_set_removed_nec 
            = vertex_handle_candidate_set_removed_nec_it->second;

        std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> intersetct_candidate;
        intersetct_candidate.reserve(
                  std::min(vertex_handle_candidate_set.size(),
                    specified_candidate_set_for_vertex.size()));
        std::set_intersection(vertex_handle_candidate_set.begin(),
                              vertex_handle_candidate_set. end (),
                       specified_candidate_set_for_vertex.begin(),
                       specified_candidate_set_for_vertex. end (),
                              std::back_inserter(intersetct_candidate));
        intersetct_candidate.swap(vertex_handle_candidate_set);
        intersetct_candidate.clear();

        intersetct_candidate.reserve(
                  std::min(vertex_handle_candidate_set_removed_nec.size(),
                               specified_candidate_set_for_vertex.size()));
        std::set_intersection(vertex_handle_candidate_set_removed_nec.begin(),
                              vertex_handle_candidate_set_removed_nec. end (),
                                  specified_candidate_set_for_vertex.begin(),
                                  specified_candidate_set_for_vertex. end (),
                              std::back_inserter(intersetct_candidate));
        intersetct_candidate.swap(vertex_handle_candidate_set_removed_nec);
        intersetct_candidate.clear();
      }
      temp_candidate_set_ptr             = std::addressof(temp_candidate_set);
      temp_candidate_set_removed_nec_ptr = std::addressof(temp_candidate_set_removed_nec);
    }
  }

  assert(temp_candidate_set_ptr);
  assert(temp_candidate_set_removed_nec_ptr);

  const CandidateSetContainer& processed_pattern_candidate_set 
                                         = *temp_candidate_set_ptr,
                               processed_pattern_candidate_set_removed_nec
                                         = *temp_candidate_set_removed_nec_ptr;

  if (restrictions.gcr()) {
    std::cout << "not implemented for gcr" << std::endl;
    return std::pair(std::vector<uint64_t>(), std::vector<float>());
  }

  // optimization for edge_literal
  bool gar_pattern_added_edge = false;
  GraphPatternType pattern_with_edge_literal(pattern);
  // there is an edge literal in literal_set, then add it as an edge
  // into the pattern for gar
  typename GUNDAM::EdgeID<GraphPatternType>::type 
    new_edge_id = GUNDAM::MaxEdgeId(pattern_with_edge_literal);
  new_edge_id++;
  for (const auto& literal : literal_set) {
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral) {
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               literal.edge_label(),
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
    if (literal.literal_type() == gar::LiteralType::kMlLiteral) {
       auto  ml_edge_label_it  = normal_to_ml_edge_label.find(literal.edge_label());
      assert(ml_edge_label_it != normal_to_ml_edge_label.end());
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               ml_edge_label_it->second,
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
  }
  gar::GraphAssociationRule<GraphPatternType,
                               DataGraphType> gar(pattern_with_edge_literal);
  gar.AddY(rhs_literal);
  // add both x and y literal into x
  for (const auto& literal : literal_set) {
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral){
      // edge literal has already been added into the
      // pattern
      continue;
    }
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    gar.AddX(literal);
  }

  auto& gar_pattern = gar.pattern();
  CandidateSetContainer gar_pattern_candidate_set;
  for (const auto& processed_pattern_candidate 
                 : processed_pattern_candidate_set) {
    auto pattern_vertex_handle = processed_pattern_candidate.first;
    // should not be nullptr
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle 
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    // should not be nullptr
    assert(gar_pattern_vertex_handle);
    auto [ gar_pattern_candidate_set_it,
           gar_pattern_candidate_set_ret ]
         = gar_pattern_candidate_set.emplace(gar_pattern_vertex_handle,
                                       processed_pattern_candidate.second);
    // should added successfully
    assert(gar_pattern_candidate_set_ret);
  }
  assert( processed_pattern_candidate_set.size()
            ==  gar_pattern_candidate_set.size() );

  CandidateSetContainer gar_pattern_candidate_set_removed_nec;
  for (const auto& processed_pattern_candidate 
                 : processed_pattern_candidate_set_removed_nec) {
    auto pattern_vertex_handle = processed_pattern_candidate.first;
    // should not be nullptr
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle 
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    // should not be nullptr
    assert(gar_pattern_vertex_handle);
    auto [ gar_pattern_candidate_set_removed_nec_it,
           gar_pattern_candidate_set_removed_nec_ret ]
         = gar_pattern_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                                   processed_pattern_candidate.second);
    // should added successfully
    assert(gar_pattern_candidate_set_removed_nec_ret);
  }
  assert( processed_pattern_candidate_set_removed_nec.size()
             == gar_pattern_candidate_set_removed_nec.size() );

  if (gar_pattern_added_edge) {
    // further refine the candidate set for more efficient processing
    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                          data_graph, 
                               gar_pattern_candidate_set_removed_nec)) {
      //return std::pair(0, _has_match::does_not_return_confidence);
    return std::pair(std::vector<uint64_t>(), std::vector<float>());
    }
    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                          data_graph, 
                                           gar_pattern_candidate_set)) {
      //return std::pair(0, _has_match::does_not_return_confidence);
    return std::pair(std::vector<uint64_t>(), std::vector<float>());
    }
  }

  std::vector<uint64_t> x_supp_counter_vec(new_literals_to_check.size(), 0);
  std::vector<uint64_t> xy_supp_counter_vec(new_literals_to_check.size(), 0);
  unsigned x_literals_should_continue = new_literals_to_check.size();
  unsigned xy_literals_should_continue = new_literals_to_check.size();


  std::function<bool(unsigned,
                        const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &x_literals_should_continue] (unsigned literal_idx,
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || x_supp_counter_vec[literal_idx] < supp_bound);
    x_supp_counter_vec[literal_idx]++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (x_supp_counter_vec[literal_idx] == supp_bound) {
      x_literals_should_continue--;
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  std::function<bool(unsigned, const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &xy_literals_should_continue](
                     unsigned literal_idx,
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    //assert(kConsiderConfidenceBound
    //    || xy_supp_counter < supp_bound);
    xy_supp_counter_vec[literal_idx]++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (xy_supp_counter_vec[literal_idx] == supp_bound) {
      xy_literals_should_continue--;
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)>
    all_satisfy_x_supp_callback 
        = [&x_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &x_literals_should_continue] (
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
   // assert(kConsiderConfidenceBound
   //     || x_supp_counter < supp_bound);
    //x_supp_counter_vec[literal_idx]++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when all reach supp_bound
      return true;
    }
    // can stop when all gars reached the support bound
    if (x_literals_should_continue == 0) {
      return false;
    }
    return true;
  };

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)> 
    all_satisfy_xy_supp_callback 
        = [&xy_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &xy_literals_should_continue](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    //assert(kConsiderConfidenceBound
    //    || xy_supp_counter < supp_bound);

    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when all reached the support bound
    if (xy_literals_should_continue == 0) {
      return false;
    }
    return true;
  };

  assert(kRhsLiteralStandAloneInfo
     == (*gar.y_literal_set()
             .begin())->stand_alone_info());

  if (!restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
    GUNDAM::Match<GraphPatternType,
                    DataGraphType> match_state;
    if (restrictions.gcr()) {
      std::cout << "not suitbale for gcr" << std::endl;
      return std::pair(std::vector<uint64_t>(), std::vector<float>());
    }
    std::vector<std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>>> match_vec_for_all;
    // gar or horn rule, use isomorphism
    auto [x_supp_vec, xy_supp_vec] = gar::GarSuppForSet<GUNDAM::MatchSemantics::kIsomorphism>(
                                          gar, data_graph, match_state,
                                          new_literals_to_check,
                                          gar_pattern_candidate_set,
                                          gar_pattern_candidate_set_removed_nec,
                                          store_match,
                                          match_vec_for_all,
                                          satisfy_x_supp_callback,
                                          satisfy_xy_supp_callback,
                                          all_satisfy_x_supp_callback,
                                          all_satisfy_xy_supp_callback,
                                          time_limit,
                                          time_limit_per_supp);
    //assert(xy_supp == xy_supp_counter);
    //assert( x_supp ==  x_supp_counter);
    if (x_supp_counter_vec.size() != new_literals_to_check.size()) {
      std::cout << "support counter is not same as literal vec" << std::endl;
    }
    for (unsigned i = 0; i < x_supp_vec.size(); i++) {
      //std::cout << "supp " << x_supp_counter_vec[i] << " supp xy " << xy_supp_counter_vec[i] << std::endl;
      if(x_supp_vec[i] != x_supp_counter_vec[i]) {
        std::cout << "not the same x_supp" << std::endl;
      }
      if(xy_supp_vec[i] != xy_supp_counter_vec[i]) {
        std::cout << "not the same xy_supp" << std::endl;
      }
    }

    if (!kConsiderConfidenceBound) {
      //assert( x_supp_counter == xy_supp_counter);
      return std::pair(xy_supp_counter_vec, std::vector<float>());
      //return std::pair(xy_supp_counter, _has_match::does_not_return_confidence);
    }
    std::vector<float> ret_conf_vec(xy_supp_counter_vec.size(), -1);
    for (unsigned i = 0; i < ret_conf_vec.size(); i++) {
      if (x_supp_counter_vec[i] == 0) {
        ret_conf_vec[i] = 0;
      } else {
        ret_conf_vec[i] = ((float) xy_supp_counter_vec[i])/ ((float)x_supp_counter_vec[i]);
      }
    }
    // if (x_supp_counter == 0) {
    //   return std::pair(xy_supp_counter, 0);
    // }
    // assert(x_supp_counter > 0);
    // return std::pair(xy_supp_counter, ((float) xy_supp_counter) 
    //                                 / ((float)  x_supp_counter));
    if (store_match) {
      match_id_vec_for_all.resize(new_literals_to_check.size());
      for (unsigned literal_idx = 0; literal_idx < match_vec_for_all.size(); literal_idx++) {
        TransferMatchToID(gar, data_graph,
                    match_vec_for_all[literal_idx], match_id_vec_for_all[literal_idx]);
      }
    }

    return std::pair(xy_supp_counter_vec, ret_conf_vec);
  }

  // specified support set 
  const auto& support_set = restrictions.support_set_for_rhs_literal(kRhsLiteralStandAloneInfo);
  assert(!support_set.empty());

  const auto kRhsLiteralInfo
         = (*gar.y_literal_set()
                .begin())->info();

  auto begin_time = std::time(NULL);

  std::vector<std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>>> match_vec_for_all;
  for (const auto& support : support_set) {
    assert(kRhsLiteralStandAloneInfo.VertexNum() == support.size());
    assert(kRhsLiteralStandAloneInfo.VertexNum() == 2);

    if (time_limit > 0 
     && time_limit < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }

    GUNDAM::Match<GraphPatternType,
                     DataGraphType> partial_match_state;

    assert(gar.pattern().FindVertex(kRhsLiteralInfo.x_id()));
    assert(gar.pattern().FindVertex(kRhsLiteralInfo.y_id()));
    partial_match_state.AddMap(gar.pattern().FindVertex(kRhsLiteralInfo.x_id()),
                               support[0]);
    partial_match_state.AddMap(gar.pattern().FindVertex(kRhsLiteralInfo.y_id()),
                               support[1]);
    
    if (restrictions.gcr()) {
      std::cout << "not suitbale for gcr" << std::endl;
      return std::pair(std::vector<uint64_t>(), std::vector<float>());;
    }
    // gar or horn rule, use isomorphism
    // gar or horn rule, use isomorphism

    // gar or horn rule, use isomorphism
    auto [x_supp_vec, xy_supp_vec] = gar::GarSuppForSet<GUNDAM::MatchSemantics::kIsomorphism>(
                                          gar, data_graph, partial_match_state,
                                          new_literals_to_check,
                                          gar_pattern_candidate_set,
                                          gar_pattern_candidate_set_removed_nec,
                                          store_match,
                                          match_vec_for_all,
                                          satisfy_x_supp_callback, 
                                          satisfy_xy_supp_callback, 
                                          all_satisfy_x_supp_callback,
                                          all_satisfy_xy_supp_callback,
                                          time_limit,
                                          time_limit_per_supp);
    //assert(xy_supp == xy_supp_counter);
    //assert( x_supp ==  x_supp_counter);
    for (unsigned i = 0; i < x_supp_vec.size(); i++) {
      if(x_supp_vec[i] != x_supp_counter_vec[i]) {
        std::cout << "not the same x_supp" << std::endl;
      }
      if(xy_supp_vec[i] != xy_supp_counter_vec[i]) {
        std::cout << "not the same xy_supp" << std::endl;
      }
    }
  }


  if (!kConsiderConfidenceBound) {
    //assert( x_supp_counter == xy_supp_counter);
    return std::pair(xy_supp_counter_vec, std::vector<float>());
    //return std::pair(xy_supp_counter, _has_match::does_not_return_confidence);
  }
  std::vector<float> ret_conf_vec(xy_supp_counter_vec.size(), -1);
  for (unsigned i = 0; i < ret_conf_vec.size(); i++) {
    if (x_supp_counter_vec[i] == 0) {
      ret_conf_vec[i] = 0;
    } else {
      ret_conf_vec[i] = ((float) xy_supp_counter_vec[i])/ ((float)x_supp_counter_vec[i]);
    }
  }

  if (store_match) {
    match_id_vec_for_all.resize(new_literals_to_check.size());
    for (unsigned literal_idx = 0; literal_idx < match_vec_for_all.size(); literal_idx++) {
      TransferMatchToID(gar, data_graph,
                  match_vec_for_all[literal_idx], match_id_vec_for_all[literal_idx]);
    }
  }

  return std::pair(xy_supp_counter_vec, ret_conf_vec);
}



template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
std::pair<std::vector<uint64_t>, std::vector<float>> HasMatchForSet(
 GraphPatternType& pattern,
 const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literal_set,
             const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal,
             const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& new_literals_to_check,
                                       DataGraphType& data_graph,
       const std::map<typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type,
                      typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type>& normal_to_ml_edge_label,
                  const CandidateSetContainer& pattern_candidate_set,
                  const CandidateSetContainer& pattern_candidate_set_removed_nec,
                  const uint64_t supp_bound,
                  const bool only_verify_is_legal,
    const GUNDAM::GraphBasicStatistics<DataGraphType>& graph_statistic,
                  const Restriction<GraphPatternType,
                                       DataGraphType>& restrictions,
                  const double time_limit,
                  const double time_limit_per_supp) {
  bool store_match = false;
  std::vector<std::vector<std::pair<std::vector<typename GUNDAM::VertexID<DataGraphType>::type>,
                        std::vector<typename GUNDAM::EdgeID<DataGraphType>::type>>>> match_id_vec_for_all;
  return HasMatchForSet(pattern,
                        literal_set,
                        rhs_literal,
                        new_literals_to_check,
                        data_graph,
                        normal_to_ml_edge_label,
                        pattern_candidate_set,
                        pattern_candidate_set_removed_nec,
                        supp_bound,
                        only_verify_is_legal,
                        graph_statistic,
                        restrictions,
                        time_limit,
                        time_limit_per_supp,
                        store_match,
                        match_id_vec_for_all);
}

template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
int64_t IncHasMatch(GraphPatternType& pattern,
 const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literal_set,
             const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal,
              DataGraphType& origin_data_graph,
    const std::vector<std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>& origin_data_graph_nec,
              DataGraphType& updated_data_graph,
    const std::vector<std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>& updated_data_graph_nec,
       const std::map<typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type,
                      typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type>& normal_to_ml_edge_label,
                  const CandidateSetContainer& pattern_origin_candidate_set,
                  const CandidateSetContainer& pattern_origin_candidate_set_removed_nec,
                  const CandidateSetContainer& pattern_updated_candidate_set,
                  const CandidateSetContainer& pattern_updated_candidate_set_removed_nec,
                  std::unordered_map<typename GUNDAM::VertexID<DataGraphType>::type,
                      std::pair<typename GUNDAM::VertexHandle<DataGraphType>::type,
                                typename GUNDAM::VertexHandle<DataGraphType>::type>>
                                id_to_origin_updated_vertex_handle,
                  std::unordered_set<typename GUNDAM::VertexID<DataGraphType>::type>
                                updated_vertex_id,
                  std::string& pivot_file_name,
                  const uint64_t supp_bound,
                  const bool only_verify_is_legal,
                  const Restriction<GraphPatternType,
                                       DataGraphType>& restrictions,
                  const double time_limit,
                  const double time_limit_per_supp) {
  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using DataGraphVertexIDType = typename GUNDAM::VertexID<DataGraphType>::type;
  const bool kConsiderConfidenceBound 
          = !only_verify_is_legal
          && restrictions.specified_confidence_bound();

  assert(literal_set.size() >= 1);
  #ifndef NDEBUG 
  bool _has_same = false;
  for (const auto& literal : literal_set) {
    if (literal == rhs_literal) {
      assert(!_has_same); // should only have one same
      _has_same = true;
    }
  }
  assert(_has_same);
  #endif // NDEBUG 

  gar::GraphAssociationRule<GraphPatternType,
                                DataGraphType> temp_gar(pattern);
  temp_gar.AddY(rhs_literal);

  const auto kRhsLiteralStandAloneInfo
          = (*temp_gar.y_literal_set()
                      .begin())->stand_alone_info();

  if (!restrictions.specified_rhs_literal_set().empty()) {
    // specified rhs literal set
    if (restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
      std::cout << "should not enter here" << std::endl;
    }
  }

  // optimization for edge_literal
  bool gar_pattern_added_edge = false;
  GraphPatternType pattern_with_edge_literal(pattern);
  // there is an edge literal in literal_set, then add it as an edge
  // into the pattern for gar
  typename GUNDAM::EdgeID<GraphPatternType>::type 
    new_edge_id = GUNDAM::MaxEdgeId(pattern_with_edge_literal);
  new_edge_id++;
  for (const auto& literal : literal_set) {
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral) {
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               literal.edge_label(),
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
    if (literal.literal_type() == gar::LiteralType::kMlLiteral) {
       auto  ml_edge_label_it  = normal_to_ml_edge_label.find(literal.edge_label());
      assert(ml_edge_label_it != normal_to_ml_edge_label.end());
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               ml_edge_label_it->second,
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
  }
  gar::GraphAssociationRule<GraphPatternType,
                               DataGraphType> gar(pattern_with_edge_literal);
  gar.AddY(rhs_literal);
  // add both x and y literal into x
  for (const auto& literal : literal_set) {
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral){
      // edge literal has already been added into the
      // pattern
      continue;
    }
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    gar.AddX(literal);
  }
  auto& gar_pattern = gar.pattern();


  CandidateSetContainer gar_pattern_origin_candidate_set;
  for (const auto& pattern_origin_candidate
                 : pattern_origin_candidate_set) {
    auto pattern_vertex_handle = pattern_origin_candidate.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_origin_candidate_set_it,
          gar_pattern_origin_candidate_set_ret]
         = gar_pattern_origin_candidate_set.emplace(gar_pattern_vertex_handle,
                                            pattern_origin_candidate.second);
    assert(gar_pattern_origin_candidate_set_ret);
  }
  assert(gar_pattern_origin_candidate_set.size()
          == pattern_origin_candidate_set.size());

  CandidateSetContainer gar_pattern_updated_candidate_set;
  for (const auto& pattern_updated_candidate
                 : pattern_updated_candidate_set) {
    auto pattern_vertex_handle = pattern_updated_candidate.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_updated_candidate_set_it,
          gar_pattern_updated_candidate_set_ret]
         = gar_pattern_updated_candidate_set.emplace(gar_pattern_vertex_handle,
                                            pattern_updated_candidate.second);
    assert(gar_pattern_updated_candidate_set_ret);
  }
  assert(gar_pattern_updated_candidate_set.size()
          == pattern_updated_candidate_set.size());

  CandidateSetContainer gar_pattern_origin_candidate_set_removed_nec;
  for (const auto& pattern_origin_candidate_removed_nec
                 : pattern_origin_candidate_set_removed_nec) {
    auto pattern_vertex_handle = pattern_origin_candidate_removed_nec.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_origin_candidate_set_removed_nec_it,
          gar_pattern_origin_candidate_set_removed_nec_ret]
         = gar_pattern_origin_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                            pattern_origin_candidate_removed_nec.second);
    assert(gar_pattern_origin_candidate_set_removed_nec_ret);
  }
  assert(gar_pattern_origin_candidate_set_removed_nec.size()
          == pattern_origin_candidate_set_removed_nec.size());

  CandidateSetContainer gar_pattern_updated_candidate_set_removed_nec;
  for (const auto& pattern_updated_candidate_removed_nec
                 : pattern_updated_candidate_set_removed_nec) {
    auto pattern_vertex_handle = pattern_updated_candidate_removed_nec.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_updated_candidate_set_removed_nec_it,
          gar_pattern_updated_candidate_set_removed_nec_ret]
         = gar_pattern_updated_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                            pattern_updated_candidate_removed_nec.second);
    assert(gar_pattern_updated_candidate_set_removed_nec_ret);
  }
  assert(gar_pattern_updated_candidate_set_removed_nec.size()
          == pattern_updated_candidate_set_removed_nec.size());


  int64_t x_supp_increase_counter = 0, x_supp_decrease_counter = 0;
  bool no_increase = false, no_decrease = false;
  if (gar_pattern_added_edge) {
    // further refine the candidate set for more efficient processing
    
    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                         origin_data_graph, 
                                            gar_pattern_origin_candidate_set)) {
      no_decrease = true;
      //std::cout << "no decrease" << std::endl;
    } else {
      gar_pattern_origin_candidate_set_removed_nec = gar_pattern_origin_candidate_set;
      GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                      DataGraphType>(
                              gar_pattern_origin_candidate_set_removed_nec,
                              origin_data_graph_nec);
    }

    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                         updated_data_graph, 
                                          gar_pattern_updated_candidate_set)) {
      no_increase = true;
      //std::cout << "no increase " << std::endl;
    } else {
      gar_pattern_updated_candidate_set_removed_nec = gar_pattern_updated_candidate_set;
      GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                      DataGraphType>(
                              gar_pattern_updated_candidate_set_removed_nec,
                              updated_data_graph_nec);
    }

    // if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
    //                                                      updated_data_graph, 
    //                            gar_pattern_updated_candidate_set_removed_nec)) {
    //   no_increase = true;
    // }
    // if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
    //                                                      updated_data_graph, 
    //                                         gar_pattern_updated_candidate_set)) {
    //   no_increase = true;
    // }
  }


  uint64_t x_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter,
           &supp_bound,
           &kConsiderConfidenceBound](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || x_supp_counter < supp_bound);
    x_supp_counter++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (x_supp_counter == supp_bound) {
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  uint64_t xy_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter,
           &supp_bound,
           &kConsiderConfidenceBound](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || xy_supp_counter < supp_bound);
    xy_supp_counter++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (xy_supp_counter == supp_bound) {
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  const auto kRhsLiteralInfo
         = (*gar.y_literal_set()
                .begin())->info();

  auto rhs_literal_x_id = kRhsLiteralInfo.x_id();
  auto rhs_literal_y_id = kRhsLiteralInfo.y_id();

  //std::cout << "start to check deleted edges" << std::endl;

  auto begin_time = std::time(NULL);

  std::unordered_map<DataGraphVertexIDType,
                      std::unordered_set<DataGraphVertexIDType>>
      deletion_affected_pivots, all_pivots;


  if (pivot_file_name.length() != 0) {
    std::ifstream pivot_file(pivot_file_name);
    if (!pivot_file.is_open()) {
      std::cout << "can not open input file " << pivot_file_name << std::endl;
      pivot_file_name.clear();
    } else {
      uint64_t data_x_pivot, data_y_pivot, affected_status;
      while (pivot_file >> data_x_pivot >> data_y_pivot >> affected_status) {
        all_pivots[data_x_pivot].insert(data_y_pivot);
        if (affected_status == 1) {
          deletion_affected_pivots[data_x_pivot].insert(data_y_pivot);
        }
      }
      pivot_file.close();
    }
  }
 
  std::map<DataGraphVertexHandleType, std::set<DataGraphVertexHandleType>>
        decrease_support_candidate_set;
  if (pivot_file_name.length() == 0) {
    // get affected pivoted nodes on origin_graph
    for (auto pattern_vertex_it = gar_pattern.VertexBegin();
            !pattern_vertex_it.IsDone();
              pattern_vertex_it++) {
      if (time_limit > 0 
        && (time_limit) < (std::time(NULL) - begin_time)) {
        // has reached time limit
        break;
      }

      for (auto pattern_out_edge_it = pattern_vertex_it->OutEdgeBegin();
              !pattern_out_edge_it.IsDone();
                pattern_out_edge_it++) {

        if (time_limit > 0 
          && (time_limit)  < (std::time(NULL) - begin_time)) {
          // has reached time limit
          break;
        }


        auto pattern_src_handle = pattern_out_edge_it->src_handle();
        auto pattern_dst_handle = pattern_out_edge_it->dst_handle();
        auto pattern_edge_label = pattern_out_edge_it->label();
        auto pattern_dst_label = pattern_dst_handle->label();

        //std::cout << "candidate set size " << gar_pattern_origin_candidate_set[pattern_src_handle].size() << std::endl;

        for (auto data_src_vertex_handle : gar_pattern_origin_candidate_set[pattern_src_handle]) {
          bool has_at_least_one_candidate = false;
          std::vector<DataGraphVertexHandleType> inc_dst_vertex_candidate;
          for (auto data_out_edge_it = data_src_vertex_handle->OutEdgeBegin(pattern_edge_label);
                  !data_out_edge_it.IsDone();
                    data_out_edge_it++) {
            auto data_dst_handle = data_out_edge_it->dst_handle();
            //std::cout << "to check an edge " << std::endl;
            const auto &attribute_handle = data_out_edge_it->FindAttribute("inc_state");
            if ((data_dst_handle->label() == pattern_dst_label)
                  && data_out_edge_it->template const_attribute<int>("inc_state") == 0) {
              //std::cout << "an edge deleted from " << data_src_vertex_handle->id() << " to " << data_dst_handle->id() << std::endl;
              has_at_least_one_candidate = true;
              inc_dst_vertex_candidate.emplace_back(data_dst_handle);
            }
          }
          if (!has_at_least_one_candidate) {
            continue;
          }
          std::vector<DataGraphVertexHandleType> inc_src_vertex_candidate;
          inc_src_vertex_candidate.emplace_back(data_src_vertex_handle);


          // CandidateSetContainer current_gar_pattern_origin_candidate_set;
          // current_gar_pattern_origin_candidate_set.emplace(pattern_src_handle, inc_src_vertex_candidate);
          // current_gar_pattern_origin_candidate_set.emplace(pattern_dst_handle, inc_dst_vertex_candidate);
          // // swap with candidate for entire graph


          // if (GUNDAM::_dp_iso::InitCandidateSetByMatchedEdges<
          //                           GUNDAM::MatchSemantics::kIsomorphism>(
          //                             gar_pattern, origin_data_graph,
          //                             current_gar_pattern_origin_candidate_set)) {
          //   if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
          //                                                 origin_data_graph, 
          //                                    current_gar_pattern_origin_candidate_set)) {
          //     GUNDAM::Match<GraphPatternType,
          //                     DataGraphType> partial_match_state;
          //     x_supp_counter = 0;
          //     xy_supp_counter = 0;
          //     auto [temp_x_supp, temp_xy_supp] = 
          //                 gar::IncGarSupp(gar, origin_data_graph, partial_match_state,
          //                 current_gar_pattern_origin_candidate_set,
          //                 decrease_support_candidate_set,
          //                 satisfy_x_supp_callback,
          //                 satisfy_xy_supp_callback,
          //                 time_limit_per_supp,
          //                 time_limit_per_supp);
          //   }
          // }

          auto current_gar_pattern_origin_candidate_set = gar_pattern_origin_candidate_set;
          // swap with candidate for entire graph
          current_gar_pattern_origin_candidate_set[pattern_src_handle].swap(inc_src_vertex_candidate);
          current_gar_pattern_origin_candidate_set[pattern_dst_handle].swap(inc_dst_vertex_candidate);


          if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                            origin_data_graph, 
                                              current_gar_pattern_origin_candidate_set)) {
              GUNDAM::Match<GraphPatternType,
                              DataGraphType> partial_match_state;
              x_supp_counter = 0;
              xy_supp_counter = 0;
              auto [temp_x_supp, temp_xy_supp] = 
                          gar::IncGarSupp(gar, origin_data_graph, partial_match_state,
                          current_gar_pattern_origin_candidate_set,
                          decrease_support_candidate_set,
                          satisfy_x_supp_callback,
                          satisfy_xy_supp_callback,
                          time_limit_per_supp,
                          time_limit_per_supp);
          }
        }
      }
    }
  } else {
    for (auto &[data_x_id, data_y_id_set] : deletion_affected_pivots) {
      auto origin_data_x_handle = id_to_origin_updated_vertex_handle[data_x_id].first;
      if (!origin_data_x_handle) {
        continue;
      }
      for (auto data_y_id : data_y_id_set) {
        auto origin_data_y_handle = id_to_origin_updated_vertex_handle[data_y_id].second;
        if (!origin_data_y_handle) {
          continue;
        }
        decrease_support_candidate_set[origin_data_x_handle].insert(origin_data_y_handle);
      }
    }
  }
  //std::cout << "######start to verify candidates#######" << std::endl;
  // for (auto &[origin_rhs_x_handle, origin_rhs_y_handle_set] : decrease_support_candidate_set) {
  //   auto data_graph_x_id = origin_rhs_x_handle->id();
  //   auto updated_graph_x_handle = id_to_origin_updated_vertex_handle[data_graph_x_id].second;
  //   if (!updated_graph_x_handle) {
  //     x_supp_decrease_counter += origin_rhs_y_handle_set.size();
  //     continue;
  //   }
  //   for (auto origin_rhs_y_handle : origin_rhs_y_handle_set) {
  //     auto data_graph_y_id = origin_rhs_y_handle->id();
  //     auto updated_graph_y_handle = id_to_origin_updated_vertex_handle[data_graph_y_id].second;
  //     if (!updated_graph_y_handle) {
  //       x_supp_decrease_counter++;
  //       continue;
  //     }
  //     decrease_updated_support_candidate_set[updated_graph_x_handle].insert(updated_graph_y_handle);
  //   }
  // }

  // auto gar_decrease_num = gar::GarSupp<GUNDAM::MatchSemantics::kIsomorphism>(
  //                                           gar, updated_data_graph, partial_match_state,
  //                                           pattern_updated_candidate_set,
  //                                           pattern_updated_candidate_set_removed_nec,
  //                                           decrease_updated_support_candidate_set,
  //                                           satisfy_x_supp_callback, 
  //                                           satisfy_xy_supp_callback, 
  //                                           time_limit_per_supp,
  //                                           time_limit_per_supp);
  // x_supp_decrease_counter += gar_decrease_num;
  
  // verify the affected pivoted node on updated graph
  begin_time = std::time(NULL);
  for (auto &[origin_rhs_x_handle, origin_rhs_y_handle_set] : decrease_support_candidate_set) {

    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
    auto data_graph_x_id = origin_rhs_x_handle->id();
    auto updated_graph_x_handle = id_to_origin_updated_vertex_handle[data_graph_x_id].second;
    if (!updated_graph_x_handle) {
  //    std::cout << "exist removed key" << std::endl;
      x_supp_decrease_counter += origin_rhs_y_handle_set.size();
      continue;
    }
    for (auto origin_rhs_y_handle : origin_rhs_y_handle_set) {
      auto data_graph_y_id = origin_rhs_y_handle->id();
      auto updated_graph_y_handle = id_to_origin_updated_vertex_handle[data_graph_y_id].second;
      if (!updated_graph_y_handle) {
//      std::cout << "exist removed  y key" << std::endl;
        x_supp_decrease_counter++;
        continue;
      }

      GUNDAM::Match<GraphPatternType,
                      DataGraphType> partial_match_state;

      // auto match_state_x_handle = gar.pattern().FindVertex(rhs_literal_x_id);
      // auto match_state_y_handle = gar.pattern().FindVertex(rhs_literal_y_id);

      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_x_id),
                                updated_graph_x_handle);
      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_y_id),
                                updated_graph_y_handle);
                              
      x_supp_counter = 0;
      xy_supp_counter = 0;


      auto [x_supp, xy_supp] = gar::GarSupp<GUNDAM::MatchSemantics::kIsomorphism>(
                                            gar, updated_data_graph, partial_match_state,
                                            gar_pattern_updated_candidate_set,
                                            gar_pattern_updated_candidate_set_removed_nec,
                                            satisfy_x_supp_callback, 
                                            satisfy_xy_supp_callback,  
                                            time_limit_per_supp,
                                            time_limit_per_supp);
      assert(xy_supp <= x_supp);
      assert(xy_supp <= 1 && x_supp <= 1);

      //std::cout << "x_supp " << x_supp << " xy supp " << xy_supp << std::endl;
      if (x_supp == 0) {
        x_supp_decrease_counter++;
      }
    }
  }
    //std::cout << "#####finish deletion verify on updated graph#####" << std::endl;
  decrease_support_candidate_set.clear();

  std::map<DataGraphVertexHandleType, std::set<DataGraphVertexHandleType>>
      increase_support_candidate_set;


  if (pivot_file_name.length() != 0) {
    for (auto &[data_x_id, data_y_id_set] : all_pivots) {
      auto updated_data_x_handle = id_to_origin_updated_vertex_handle[data_x_id].second;
      if (!updated_data_x_handle) {
        continue;
      }
      for (auto data_y_id : data_y_id_set) {
        auto updated_data_y_handle = id_to_origin_updated_vertex_handle[data_y_id].second;
        if (!updated_data_y_handle) {
          continue;
        }
        increase_support_candidate_set[updated_data_x_handle].insert(updated_data_y_handle);
      }
    }
  }

// get affected pivoted nodes on updated graph
  begin_time = std::time(NULL);
  for (auto pattern_vertex_it = gar_pattern.VertexBegin();
           !pattern_vertex_it.IsDone();
            pattern_vertex_it++) {
    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
    for (auto pattern_out_edge_it = pattern_vertex_it->OutEdgeBegin();
             !pattern_out_edge_it.IsDone();
              pattern_out_edge_it++) {
    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
      auto pattern_src_handle = pattern_out_edge_it->src_handle();
      auto pattern_dst_handle = pattern_out_edge_it->dst_handle();
      auto pattern_edge_label = pattern_out_edge_it->label();
      auto pattern_dst_label = pattern_dst_handle->label();

      for (auto data_src_vertex_handle : gar_pattern_updated_candidate_set[pattern_src_handle]) {
        bool has_at_least_one_candidate = false;
        std::vector<DataGraphVertexHandleType> inc_dst_vertex_candidate;
        for (auto data_out_edge_it = data_src_vertex_handle->OutEdgeBegin(pattern_edge_label);
                 !data_out_edge_it.IsDone();
                  data_out_edge_it++) {
          auto data_dst_handle = data_out_edge_it->dst_handle();
          if ((data_dst_handle->label() == pattern_dst_label)
                && data_out_edge_it->template const_attribute<int>("inc_state") == 2) {
            has_at_least_one_candidate = true;
            inc_dst_vertex_candidate.emplace_back(data_dst_handle);
          }
        }
        if (!has_at_least_one_candidate) {
          continue;
        }
        std::vector<DataGraphVertexHandleType> inc_src_vertex_candidate;
        inc_src_vertex_candidate.emplace_back(data_src_vertex_handle);


        // CandidateSetContainer current_gar_pattern_updated_candidate_set;
        // current_gar_pattern_updated_candidate_set.emplace(pattern_src_handle, inc_src_vertex_candidate);
        // current_gar_pattern_updated_candidate_set.emplace(pattern_dst_handle, inc_dst_vertex_candidate);
        // swap with candidate for entire graph


        // if (GUNDAM::_dp_iso::InitCandidateSetByMatchedEdges<
        //                           GUNDAM::MatchSemantics::kIsomorphism>(
        //                             gar_pattern, updated_data_graph,
        //                             current_gar_pattern_updated_candidate_set)) {
        //   if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
        //                                                 updated_data_graph, 
        //                                    current_gar_pattern_updated_candidate_set)) {
        //     GUNDAM::Match<GraphPatternType,
        //                     DataGraphType> partial_match_state;
        //     x_supp_counter = 0;
        //     xy_supp_counter = 0;
        //     auto [temp_x_supp, temp_xy_supp] = 
        //                 gar::IncGarSupp(gar, updated_data_graph, partial_match_state,
        //                 current_gar_pattern_updated_candidate_set,
        //                 increase_support_candidate_set,
        //                 satisfy_x_supp_callback,
        //                 satisfy_xy_supp_callback,
        //                 time_limit_per_supp,
        //                 time_limit_per_supp);
        //   }
        // }

        auto current_gar_pattern_updated_candidate_set = gar_pattern_updated_candidate_set;
        current_gar_pattern_updated_candidate_set[pattern_src_handle].swap(inc_src_vertex_candidate);
        current_gar_pattern_updated_candidate_set[pattern_dst_handle].swap(inc_dst_vertex_candidate);


       if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                        updated_data_graph, 
                                           current_gar_pattern_updated_candidate_set)) {
          GUNDAM::Match<GraphPatternType,
                          DataGraphType> partial_match_state;
          x_supp_counter = 0;
          xy_supp_counter = 0;
          auto [temp_x_supp, temp_xy_supp] = 
                      gar::IncGarSupp(gar, updated_data_graph, partial_match_state,
                      current_gar_pattern_updated_candidate_set,
                      increase_support_candidate_set,
                      satisfy_x_supp_callback,
                      satisfy_xy_supp_callback,
                      time_limit_per_supp,
                      time_limit_per_supp);
       }
      }
    }
  }



// verify the affected pivoted node on origin graph
  begin_time = std::time(NULL);
  for (auto &[updated_rhs_x_handle, updated_rhs_y_handle_set] : increase_support_candidate_set) {
    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
    auto data_graph_x_id = updated_rhs_x_handle->id();
    auto origin_graph_x_handle = id_to_origin_updated_vertex_handle[data_graph_x_id].first;
    if (!origin_graph_x_handle) {
      x_supp_increase_counter += updated_rhs_y_handle_set.size();
      continue;
    }
    for (auto updated_rhs_y_handle : updated_rhs_y_handle_set) {
      auto data_graph_y_id = updated_rhs_y_handle->id();
      auto origin_graph_y_handle = id_to_origin_updated_vertex_handle[data_graph_y_id].first;
      if (!origin_graph_y_handle) {
        x_supp_increase_counter++;
        continue;
      }
      if((all_pivots.find(data_graph_x_id) != all_pivots.end())
          && (all_pivots[data_graph_x_id].find(data_graph_y_id)
              != all_pivots[data_graph_x_id].end())) {
        continue;
      }

      GUNDAM::Match<GraphPatternType,
                      DataGraphType> partial_match_state;

      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_x_id),
                                origin_graph_x_handle);
      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_y_id),
                                origin_graph_y_handle);
      x_supp_counter = 0;
      xy_supp_counter = 0;
      auto [x_supp, xy_supp] = gar::GarSupp<GUNDAM::MatchSemantics::kIsomorphism>(
                                            gar, origin_data_graph, partial_match_state,
                                            gar_pattern_origin_candidate_set,
                                            gar_pattern_origin_candidate_set_removed_nec,
                                            satisfy_x_supp_callback, 
                                            satisfy_xy_supp_callback, 
                                            time_limit_per_supp,
                                            time_limit_per_supp);
      assert(xy_supp <= x_supp);
      assert(xy_supp <= 1 && x_supp <= 1);

      
      if (x_supp == 0) {
        x_supp_increase_counter++;
      }
    }
  }

  return (int64_t)x_supp_increase_counter - (int64_t)x_supp_decrease_counter;
}


template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
std::vector<int64_t> IncHasMatchForSet(GraphPatternType& pattern,
 const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literal_set,
             const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal,
             const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literals_to_check,
              DataGraphType& origin_data_graph,
    const std::vector<std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>& origin_data_graph_nec,
              DataGraphType& updated_data_graph,
    const std::vector<std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>& updated_data_graph_nec,
       const std::map<typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type,
                      typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type>& normal_to_ml_edge_label,
                  const CandidateSetContainer& pattern_origin_candidate_set,
                  const CandidateSetContainer& pattern_origin_candidate_set_removed_nec,
                  const CandidateSetContainer& pattern_updated_candidate_set,
                  const CandidateSetContainer& pattern_updated_candidate_set_removed_nec,
                  std::unordered_map<typename GUNDAM::VertexID<DataGraphType>::type,
                      std::pair<typename GUNDAM::VertexHandle<DataGraphType>::type,
                                typename GUNDAM::VertexHandle<DataGraphType>::type>>
                                id_to_origin_updated_vertex_handle,
                  std::unordered_set<typename GUNDAM::VertexID<DataGraphType>::type>
                                updated_vertex_id,
                  std::vector<std::string>& pivot_file_name_vec,
                  const uint64_t supp_bound,
                  const bool only_verify_is_legal,
                  const Restriction<GraphPatternType,
                                       DataGraphType>& restrictions,
                  const double time_limit,
                  const double time_limit_per_supp) {
  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using DataGraphVertexIDType = typename GUNDAM::VertexID<DataGraphType>::type;
  const bool kConsiderConfidenceBound 
          = !only_verify_is_legal
          && restrictions.specified_confidence_bound();

  assert(literal_set.size() >= 1);
  #ifndef NDEBUG 
  bool _has_same = false;
  for (const auto& literal : literal_set) {
    if (literal == rhs_literal) {
      assert(!_has_same); // should only have one same
      _has_same = true;
    }
  }
  assert(_has_same);
  #endif // NDEBUG 

  gar::GraphAssociationRule<GraphPatternType,
                                DataGraphType> temp_gar(pattern);
  temp_gar.AddY(rhs_literal);

  const auto kRhsLiteralStandAloneInfo
          = (*temp_gar.y_literal_set()
                      .begin())->stand_alone_info();

  if (!restrictions.specified_rhs_literal_set().empty()) {
    // specified rhs literal set
    if (restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
      std::cout << "should not enter here" << std::endl;
    }
  }

  // optimization for edge_literal
  bool gar_pattern_added_edge = false;
  GraphPatternType pattern_with_edge_literal(pattern);
  // there is an edge literal in literal_set, then add it as an edge
  // into the pattern for gar
  typename GUNDAM::EdgeID<GraphPatternType>::type 
    new_edge_id = GUNDAM::MaxEdgeId(pattern_with_edge_literal);
  new_edge_id++;
  for (const auto& literal : literal_set) {
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral) {
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               literal.edge_label(),
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
    if (literal.literal_type() == gar::LiteralType::kMlLiteral) {
       auto  ml_edge_label_it  = normal_to_ml_edge_label.find(literal.edge_label());
      assert(ml_edge_label_it != normal_to_ml_edge_label.end());
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               ml_edge_label_it->second,
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
  }
  gar::GraphAssociationRule<GraphPatternType,
                               DataGraphType> gar(pattern_with_edge_literal);
  gar.AddY(rhs_literal);
  // add both x and y literal into x
  for (const auto& literal : literal_set) {
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral){
      // edge literal has already been added into the
      // pattern
      continue;
    }
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    gar.AddX(literal);
  }
  auto& gar_pattern = gar.pattern();


  CandidateSetContainer gar_pattern_origin_candidate_set;
  for (const auto& pattern_origin_candidate
                 : pattern_origin_candidate_set) {
    auto pattern_vertex_handle = pattern_origin_candidate.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_origin_candidate_set_it,
          gar_pattern_origin_candidate_set_ret]
         = gar_pattern_origin_candidate_set.emplace(gar_pattern_vertex_handle,
                                            pattern_origin_candidate.second);
    assert(gar_pattern_origin_candidate_set_ret);
  }
  assert(gar_pattern_origin_candidate_set.size()
          == pattern_origin_candidate_set.size());

  CandidateSetContainer gar_pattern_updated_candidate_set;
  for (const auto& pattern_updated_candidate
                 : pattern_updated_candidate_set) {
    auto pattern_vertex_handle = pattern_updated_candidate.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_updated_candidate_set_it,
          gar_pattern_updated_candidate_set_ret]
         = gar_pattern_updated_candidate_set.emplace(gar_pattern_vertex_handle,
                                            pattern_updated_candidate.second);
    assert(gar_pattern_updated_candidate_set_ret);
  }
  assert(gar_pattern_updated_candidate_set.size()
          == pattern_updated_candidate_set.size());

  CandidateSetContainer gar_pattern_origin_candidate_set_removed_nec;
  for (const auto& pattern_origin_candidate_removed_nec
                 : pattern_origin_candidate_set_removed_nec) {
    auto pattern_vertex_handle = pattern_origin_candidate_removed_nec.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_origin_candidate_set_removed_nec_it,
          gar_pattern_origin_candidate_set_removed_nec_ret]
         = gar_pattern_origin_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                            pattern_origin_candidate_removed_nec.second);
    assert(gar_pattern_origin_candidate_set_removed_nec_ret);
  }
  assert(gar_pattern_origin_candidate_set_removed_nec.size()
          == pattern_origin_candidate_set_removed_nec.size());

  CandidateSetContainer gar_pattern_updated_candidate_set_removed_nec;
  for (const auto& pattern_updated_candidate_removed_nec
                 : pattern_updated_candidate_set_removed_nec) {
    auto pattern_vertex_handle = pattern_updated_candidate_removed_nec.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_updated_candidate_set_removed_nec_it,
          gar_pattern_updated_candidate_set_removed_nec_ret]
         = gar_pattern_updated_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                            pattern_updated_candidate_removed_nec.second);
    assert(gar_pattern_updated_candidate_set_removed_nec_ret);
  }
  assert(gar_pattern_updated_candidate_set_removed_nec.size()
          == pattern_updated_candidate_set_removed_nec.size());


  std::vector<int64_t> x_supp_increase_counter(literals_to_check.size(), 0),
                       x_supp_decrease_counter(literals_to_check.size(), 0);
  bool no_increase = false, no_decrease = false;
  if (gar_pattern_added_edge) {
    // further refine the candidate set for more efficient processing
    
    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                         origin_data_graph, 
                                            gar_pattern_origin_candidate_set)) {
      no_decrease = true;
      //std::cout << "no decrease" << std::endl;
    } else {
      gar_pattern_origin_candidate_set_removed_nec = gar_pattern_origin_candidate_set;
      GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                      DataGraphType>(
                              gar_pattern_origin_candidate_set_removed_nec,
                              origin_data_graph_nec);
    }

    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                         updated_data_graph, 
                                          gar_pattern_updated_candidate_set)) {
      no_increase = true;
      //std::cout << "no increase " << std::endl;
    } else {
      gar_pattern_updated_candidate_set_removed_nec = gar_pattern_updated_candidate_set;
      GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                      DataGraphType>(
                              gar_pattern_updated_candidate_set_removed_nec,
                              updated_data_graph_nec);
    }

    // if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
    //                                                      updated_data_graph, 
    //                            gar_pattern_updated_candidate_set_removed_nec)) {
    //   no_increase = true;
    // }
    // if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
    //                                                      updated_data_graph, 
    //                                         gar_pattern_updated_candidate_set)) {
    //   no_increase = true;
    // }
  }


  std::vector<uint64_t> x_supp_counter_vec(literals_to_check.size(), 0);
  std::vector<uint64_t> xy_supp_counter_vec(literals_to_check.size(), 0);
  unsigned x_literals_should_continue = literals_to_check.size();
  unsigned xy_literals_should_continue = literals_to_check.size();

  std::vector<uint64_t> global_x_supp_counter_vec(literals_to_check.size(), 0);
  std::vector<uint64_t> global_xy_supp_counter_vec(literals_to_check.size(), 0);
  unsigned global_x_literals_should_continue = literals_to_check.size();
  unsigned global_xy_literals_should_continue = literals_to_check.size();


  std::function<bool(unsigned,
                        const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &x_literals_should_continue] (unsigned literal_idx,
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || x_supp_counter_vec[literal_idx] < supp_bound);
    x_supp_counter_vec[literal_idx]++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (x_supp_counter_vec[literal_idx] == supp_bound) {
      x_literals_should_continue--;
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  std::function<bool(unsigned, const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &xy_literals_should_continue](
                     unsigned literal_idx,
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    //assert(kConsiderConfidenceBound
    //    || xy_supp_counter < supp_bound);
    xy_supp_counter_vec[literal_idx]++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (xy_supp_counter_vec[literal_idx] == supp_bound) {
      xy_literals_should_continue--;
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)>
    all_satisfy_x_supp_callback 
        = [&x_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &x_literals_should_continue] (
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
   // assert(kConsiderConfidenceBound
   //     || x_supp_counter < supp_bound);
    //x_supp_counter_vec[literal_idx]++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when all reach supp_bound
      return true;
    }
    // can stop when all gars reached the support bound
    if (x_literals_should_continue == 0) {
      return false;
    }
    return true;
  };

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)> 
    all_satisfy_xy_supp_callback 
        = [&xy_supp_counter_vec,
           &supp_bound,
           &kConsiderConfidenceBound,
           &xy_literals_should_continue](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    //assert(kConsiderConfidenceBound
    //    || xy_supp_counter < supp_bound);

    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when all reached the support bound
    if (xy_literals_should_continue == 0) {
      return false;
    }
    return true;
  };



  const auto kRhsLiteralInfo
         = (*gar.y_literal_set()
                .begin())->info();

  auto rhs_literal_x_id = kRhsLiteralInfo.x_id();
  auto rhs_literal_y_id = kRhsLiteralInfo.y_id();

  //std::cout << "start to check deleted edges" << std::endl;

  auto begin_time = std::time(NULL);

  std::unordered_map<DataGraphVertexIDType,
                      std::unordered_map<DataGraphVertexIDType, std::set<unsigned>>>
      deletion_affected_pivots, all_pivots;

  for (unsigned i = 0; i < literals_to_check.size(); i++) {
    std::string pivot_file_name = pivot_file_name_vec[i];
    std::ifstream pivot_file(pivot_file_name);
    if (!pivot_file.is_open()) {
      std::cout << "can not open input file " << pivot_file_name << std::endl;
      pivot_file_name.clear();
    } else {
      uint64_t data_x_pivot, data_y_pivot, affected_status;
      while (pivot_file >> data_x_pivot >> data_y_pivot >> affected_status) {
        all_pivots[data_x_pivot][data_y_pivot].insert(i);
        if (affected_status == 1) {
          deletion_affected_pivots[data_x_pivot][data_y_pivot].insert(i);
        }
      }
      pivot_file.close();
    }
  }

 
  std::map<DataGraphVertexHandleType, std::map<DataGraphVertexHandleType, std::set<unsigned>>>
        decrease_support_candidate_set;

    for (auto &[data_x_id, data_y_id_map] : deletion_affected_pivots) {
      auto origin_data_x_handle = id_to_origin_updated_vertex_handle[data_x_id].first;
      if (!origin_data_x_handle) {
        continue;
      }
      for (auto& [data_y_id, related_literals] : data_y_id_map) {
        auto origin_data_y_handle = id_to_origin_updated_vertex_handle[data_y_id].second;
        if (!origin_data_y_handle) {
          continue;
        }
        decrease_support_candidate_set[origin_data_x_handle][origin_data_y_handle] = related_literals;
      }
    }
  
  // verify the affected pivoted node on updated graph
  begin_time = std::time(NULL);
  for (auto &[origin_rhs_x_handle, origin_rhs_y_handle_map] : decrease_support_candidate_set) {
    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
    auto data_graph_x_id = origin_rhs_x_handle->id();
    auto updated_graph_x_handle = id_to_origin_updated_vertex_handle[data_graph_x_id].second;
    if (!updated_graph_x_handle) {
  //    std::cout << "exist removed key" << std::endl;
      for (auto &[y_handle, ids] : origin_rhs_y_handle_map) {
        for (auto &id : ids) {
          x_supp_decrease_counter[id]++;
        }
      }
      continue;
    }
    for (auto &[origin_rhs_y_handle, ids] : origin_rhs_y_handle_map) {
      auto data_graph_y_id = origin_rhs_y_handle->id();
      auto updated_graph_y_handle = id_to_origin_updated_vertex_handle[data_graph_y_id].second;
      if (!updated_graph_y_handle) {
//      std::cout << "exist removed  y key" << std::endl;
        for (auto &id : ids) {
          x_supp_decrease_counter[id]++;
        }
        continue;
      }

      std::map<int, int> local_literal_id2global;
      std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>> temp_literal_to_check;
      for (auto id : ids) {
        temp_literal_to_check.emplace_back(literals_to_check[id]);
        local_literal_id2global[temp_literal_to_check.size() - 1] = id;
      }

      GUNDAM::Match<GraphPatternType,
                      DataGraphType> partial_match_state;

      // auto match_state_x_handle = gar.pattern().FindVertex(rhs_literal_x_id);
      // auto match_state_y_handle = gar.pattern().FindVertex(rhs_literal_y_id);

      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_x_id),
                                updated_graph_x_handle);
      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_y_id),
                                updated_graph_y_handle);


      x_supp_counter_vec.clear();
      x_supp_counter_vec.resize(temp_literal_to_check.size(), 0);
      xy_supp_counter_vec.clear();
      xy_supp_counter_vec.resize(temp_literal_to_check.size(), 0);
      x_literals_should_continue = temp_literal_to_check.size();
      xy_literals_should_continue = temp_literal_to_check.size();



      auto [x_supp, xy_supp] = gar::GarSuppForSet<GUNDAM::MatchSemantics::kIsomorphism>(
                                            gar, updated_data_graph, partial_match_state,
                                            temp_literal_to_check,
                                            gar_pattern_updated_candidate_set,
                                            gar_pattern_updated_candidate_set_removed_nec,
                                            satisfy_x_supp_callback,
                                            satisfy_xy_supp_callback,
                                            all_satisfy_x_supp_callback,
                                            all_satisfy_xy_supp_callback,
                                            time_limit_per_supp,
                                            time_limit_per_supp);

      for (unsigned i = 0; i < xy_supp.size(); i++) {
        if (xy_supp[i] == 0) {
          auto global_id = local_literal_id2global[i];
          x_supp_decrease_counter[global_id]++;
        }
      }

      //std::cout << "x_supp " << x_supp << " xy supp " << xy_supp << std::endl;
      // if (x_supp == 0) {
      //   x_supp_decrease_counter++;
      // }
    }
  }
    //std::cout << "#####finish deletion verify on updated graph#####" << std::endl;
  decrease_support_candidate_set.clear();

  std::map<DataGraphVertexHandleType, std::map<DataGraphVertexHandleType, std::set<unsigned>>>
        increase_support_candidate_set;



  for (unsigned i = 0; i < literals_to_check.size(); i++) {
    for (auto &[data_x_id, data_y_id_map] : all_pivots) {
      auto updated_data_x_handle = id_to_origin_updated_vertex_handle[data_x_id].second;
      if (!updated_data_x_handle) {
        continue;
      }
      for (auto &[data_y_id, ids]: data_y_id_map) {
        auto updated_data_y_handle = id_to_origin_updated_vertex_handle[data_y_id].second;
        if (!updated_data_y_handle) {
          continue;
        }
        for (auto &id : ids) {
            increase_support_candidate_set[updated_data_x_handle][updated_data_y_handle].insert(id);
        }

      }
    }
  }


// get affected pivoted nodes on updated graph
  begin_time = std::time(NULL);
  for (auto pattern_vertex_it = gar_pattern.VertexBegin();
           !pattern_vertex_it.IsDone();
            pattern_vertex_it++) {
    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
    for (auto pattern_out_edge_it = pattern_vertex_it->OutEdgeBegin();
             !pattern_out_edge_it.IsDone();
              pattern_out_edge_it++) {
    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
      auto pattern_src_handle = pattern_out_edge_it->src_handle();
      auto pattern_dst_handle = pattern_out_edge_it->dst_handle();
      auto pattern_edge_label = pattern_out_edge_it->label();
      auto pattern_dst_label = pattern_dst_handle->label();

      for (auto data_src_vertex_handle : gar_pattern_updated_candidate_set[pattern_src_handle]) {
        bool has_at_least_one_candidate = false;
        std::vector<DataGraphVertexHandleType> inc_dst_vertex_candidate;
        for (auto data_out_edge_it = data_src_vertex_handle->OutEdgeBegin(pattern_edge_label);
                 !data_out_edge_it.IsDone();
                  data_out_edge_it++) {
          auto data_dst_handle = data_out_edge_it->dst_handle();
          if ((data_dst_handle->label() == pattern_dst_label)
                && data_out_edge_it->template const_attribute<int>("inc_state") == 2) {
            has_at_least_one_candidate = true;
            inc_dst_vertex_candidate.emplace_back(data_dst_handle);
          }
        }
        if (!has_at_least_one_candidate) {
          continue;
        }
        std::vector<DataGraphVertexHandleType> inc_src_vertex_candidate;
        inc_src_vertex_candidate.emplace_back(data_src_vertex_handle);

        auto current_gar_pattern_updated_candidate_set = gar_pattern_updated_candidate_set;
        current_gar_pattern_updated_candidate_set[pattern_src_handle].swap(inc_src_vertex_candidate);
        current_gar_pattern_updated_candidate_set[pattern_dst_handle].swap(inc_dst_vertex_candidate);


        x_supp_counter_vec.clear();
        x_supp_counter_vec.resize(literals_to_check.size(), 0);
        xy_supp_counter_vec.clear();
        xy_supp_counter_vec.resize(literals_to_check.size(), 0);
        x_literals_should_continue = literals_to_check.size();
        xy_literals_should_continue = literals_to_check.size();


       if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                        updated_data_graph, 
                                           current_gar_pattern_updated_candidate_set)) {
          GUNDAM::Match<GraphPatternType,
                          DataGraphType> partial_match_state;
          // x_supp_counter = 0;
          // xy_supp_counter = 0;
          auto [temp_x_supp_vec, temp_xy_supp_vec] = 
                      gar::IncGarSuppForSet<GUNDAM::MatchSemantics::kIsomorphism>(
                      gar,
                      updated_data_graph, partial_match_state,
                      literals_to_check,
                      increase_support_candidate_set,                      
                      current_gar_pattern_updated_candidate_set,                      
                      current_gar_pattern_updated_candidate_set,
                      satisfy_x_supp_callback,
                      satisfy_xy_supp_callback,
                      all_satisfy_x_supp_callback,
                      all_satisfy_xy_supp_callback,                      
                      time_limit,
                      time_limit_per_supp);
        }
      }
    }
  }



// verify the affected pivoted node on origin graph
  begin_time = std::time(NULL);
  for (auto &[updated_rhs_x_handle, updated_rhs_y_handle_map] : increase_support_candidate_set) {
    if (time_limit > 0 
       && (time_limit)  < (std::time(NULL) - begin_time)) {
      // has reached time limit
      break;
    }
    auto data_graph_x_id = updated_rhs_x_handle->id();
    auto origin_graph_x_handle = id_to_origin_updated_vertex_handle[data_graph_x_id].first;
    if (!origin_graph_x_handle) {
      for (auto &[rhs_y_handle, ids] : updated_rhs_y_handle_map) {
        for (auto id : ids) {
          x_supp_increase_counter[id]++;
        }
      }
      continue;
    }
    for (auto &[updated_rhs_y_handle, ids] : updated_rhs_y_handle_map) {
      auto data_graph_y_id = updated_rhs_y_handle->id();
      auto origin_graph_y_handle = id_to_origin_updated_vertex_handle[data_graph_y_id].first;
      if (!origin_graph_y_handle) {
        for (auto &id : ids) {
          x_supp_increase_counter[id]++;
        }
        continue;
      }
      if((all_pivots.find(data_graph_x_id) != all_pivots.end())
          && (all_pivots[data_graph_x_id].find(data_graph_y_id)
              != all_pivots[data_graph_x_id].end())
          && (all_pivots[data_graph_x_id][data_graph_y_id]
              == ids)) {
        continue;
      }
      std::set<unsigned> ids_to_check;
      if ((all_pivots.find(data_graph_x_id) != all_pivots.end())
          && (all_pivots[data_graph_x_id].find(data_graph_y_id)
              != all_pivots[data_graph_x_id].end())) {
        auto& exist_ids = all_pivots[data_graph_x_id][data_graph_y_id];
        for (auto id : ids) {
          if (exist_ids.find(id) == exist_ids.end()) {
            ids_to_check.insert(id);
          }
        }
      } else {
        ids_to_check = ids;
      }

      std::map<int, int> local_literal_id2global;
      std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>> temp_literal_to_check;
      for (auto id : ids_to_check) {
        temp_literal_to_check.emplace_back(literals_to_check[id]);
        local_literal_id2global[temp_literal_to_check.size() - 1] = id;
      }

      x_supp_counter_vec.clear();
      x_supp_counter_vec.resize(temp_literal_to_check.size(), 0);
      xy_supp_counter_vec.clear();
      xy_supp_counter_vec.resize(temp_literal_to_check.size(), 0);
      x_literals_should_continue = temp_literal_to_check.size();
      xy_literals_should_continue = temp_literal_to_check.size();


      GUNDAM::Match<GraphPatternType,
                      DataGraphType> partial_match_state;

      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_x_id),
                                origin_graph_x_handle);
      partial_match_state.AddMap(gar.pattern().FindVertex(rhs_literal_y_id),
                                origin_graph_y_handle);


      auto [x_supp, xy_supp] = gar::GarSuppForSet<GUNDAM::MatchSemantics::kIsomorphism>(
                                            gar, origin_data_graph, partial_match_state,
                                            temp_literal_to_check,
                                            gar_pattern_origin_candidate_set,
                                            gar_pattern_origin_candidate_set_removed_nec,
                                            satisfy_x_supp_callback, 
                                            satisfy_xy_supp_callback,
                                            all_satisfy_x_supp_callback,
                                            all_satisfy_xy_supp_callback,
                                            time_limit_per_supp,
                                            time_limit_per_supp);
      assert(xy_supp <= x_supp);
      assert(xy_supp <= 1 && x_supp <= 1);

      
      for (unsigned i = 0; i < xy_supp.size(); i++) {
        if (xy_supp[i] == 0) {
          auto global_id = local_literal_id2global[i];
          x_supp_increase_counter[global_id]++;
        }
      }
    }
  }

  std::vector<int64_t> varied_supp(x_supp_increase_counter.size(), 0);
  for (unsigned i = 0; i < x_supp_increase_counter.size(); i++) {
    varied_supp[i] = (int64_t)x_supp_increase_counter[i] - (int64_t)x_supp_decrease_counter[i];
  }

  return varied_supp;
}

template <typename GraphPatternType,
          typename    DataGraphType, 
          typename CandidateSetContainer>
bool IncHasMatchByEstimate(GraphPatternType& pattern,
 const std::vector<gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>>& literal_set,
             const gar::LiteralInfo<GraphPatternType,
                                       DataGraphType>&  rhs_literal,
              DataGraphType& origin_data_graph,
    const std::vector<std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>& origin_data_graph_nec,
              DataGraphType& updated_data_graph,
    const std::vector<std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>& updated_data_graph_nec,
       const std::map<typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type,
                      typename GUNDAM::EdgeLabel<
                                       DataGraphType>::type>& normal_to_ml_edge_label,
                  const CandidateSetContainer& pattern_origin_candidate_set,
                  const CandidateSetContainer& pattern_origin_candidate_set_removed_nec,
                  const CandidateSetContainer& pattern_updated_candidate_set,
                  const CandidateSetContainer& pattern_updated_candidate_set_removed_nec,
                  std::unordered_map<typename GUNDAM::VertexID<DataGraphType>::type,
                      std::pair<typename GUNDAM::VertexHandle<DataGraphType>::type,
                                typename GUNDAM::VertexHandle<DataGraphType>::type>>
                                id_to_origin_updated_vertex_handle,
                  std::unordered_set<typename GUNDAM::VertexID<DataGraphType>::type>
                                updated_vertex_id,
                  std::string& pivot_file_name,
                  uint64_t support,
                  const uint64_t supp_bound,
                  double probability,
                  double kThreshold,
                  int64_t &estimated_support,
                  const bool only_verify_is_legal,
                  const Restriction<GraphPatternType,
                                       DataGraphType>& restrictions,
                  const double time_limit,
                  const double time_limit_per_supp) {
  using GraphPatternVertexHandleType = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandleType = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using DataGraphVertexIDType = typename GUNDAM::VertexID<DataGraphType>::type;
  const bool kConsiderConfidenceBound 
          = !only_verify_is_legal
          && restrictions.specified_confidence_bound();

  assert(literal_set.size() >= 1);
  #ifndef NDEBUG 
  bool _has_same = false;
  for (const auto& literal : literal_set) {
    if (literal == rhs_literal) {
      assert(!_has_same); // should only have one same
      _has_same = true;
    }
  }
  assert(_has_same);
  #endif // NDEBUG 

  gar::GraphAssociationRule<GraphPatternType,
                                DataGraphType> temp_gar(pattern);
  temp_gar.AddY(rhs_literal);

  const auto kRhsLiteralStandAloneInfo
          = (*temp_gar.y_literal_set()
                      .begin())->stand_alone_info();


  if (!restrictions.specified_rhs_literal_set().empty()) {
    // specified rhs literal set
    if (restrictions.specified_support_set_for_rhs_literal(kRhsLiteralStandAloneInfo)) {
      std::cout << "should not enter here" << std::endl;
    }
  }

  // optimization for edge_literal
  bool gar_pattern_added_edge = false;
  GraphPatternType pattern_with_edge_literal(pattern);
  // there is an edge literal in literal_set, then add it as an edge
  // into the pattern for gar
  typename GUNDAM::EdgeID<GraphPatternType>::type 
    new_edge_id = GUNDAM::MaxEdgeId(pattern_with_edge_literal);
  new_edge_id++;
  for (const auto& literal : literal_set) {
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral) {
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               literal.edge_label(),
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
    if (literal.literal_type() == gar::LiteralType::kMlLiteral) {
       auto  ml_edge_label_it  = normal_to_ml_edge_label.find(literal.edge_label());
      assert(ml_edge_label_it != normal_to_ml_edge_label.end());
      auto [ edge_handle, edge_ret ]
           = pattern_with_edge_literal.AddEdge(literal.x_id(),
                                               literal.y_id(),
                                               ml_edge_label_it->second,
                                               new_edge_id++);
      assert(edge_ret);
      gar_pattern_added_edge = true;
    }
  }
  gar::GraphAssociationRule<GraphPatternType,
                               DataGraphType> gar(pattern_with_edge_literal);
  gar.AddY(rhs_literal);
  // add both x and y literal into x
  for (const auto& literal : literal_set) {
    if (literal.literal_type() == gar::LiteralType::kEdgeLiteral){
      // edge literal has already been added into the
      // pattern
      continue;
    }
    if (kConsiderConfidenceBound) {
      // does not add rhs literal to lhs
      if (literal == rhs_literal) {
        continue;
      }
    }
    gar.AddX(literal);
  }
  auto& gar_pattern = gar.pattern();


  CandidateSetContainer gar_pattern_origin_candidate_set;
  for (const auto& pattern_origin_candidate
                 : pattern_origin_candidate_set) {
    auto pattern_vertex_handle = pattern_origin_candidate.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_origin_candidate_set_it,
          gar_pattern_origin_candidate_set_ret]
         = gar_pattern_origin_candidate_set.emplace(gar_pattern_vertex_handle,
                                            pattern_origin_candidate.second);
    assert(gar_pattern_origin_candidate_set_ret);
  }
  assert(gar_pattern_origin_candidate_set.size()
          == pattern_origin_candidate_set.size());

  CandidateSetContainer gar_pattern_updated_candidate_set;
  for (const auto& pattern_updated_candidate
                 : pattern_updated_candidate_set) {
    auto pattern_vertex_handle = pattern_updated_candidate.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_updated_candidate_set_it,
          gar_pattern_updated_candidate_set_ret]
         = gar_pattern_updated_candidate_set.emplace(gar_pattern_vertex_handle,
                                            pattern_updated_candidate.second);
    assert(gar_pattern_updated_candidate_set_ret);
  }
  assert(gar_pattern_updated_candidate_set.size()
          == pattern_updated_candidate_set.size());

  CandidateSetContainer gar_pattern_origin_candidate_set_removed_nec;
  for (const auto& pattern_origin_candidate_removed_nec
                 : pattern_origin_candidate_set_removed_nec) {
    auto pattern_vertex_handle = pattern_origin_candidate_removed_nec.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_origin_candidate_set_removed_nec_it,
          gar_pattern_origin_candidate_set_removed_nec_ret]
         = gar_pattern_origin_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                            pattern_origin_candidate_removed_nec.second);
    assert(gar_pattern_origin_candidate_set_removed_nec_ret);
  }
  assert(gar_pattern_origin_candidate_set_removed_nec.size()
          == pattern_origin_candidate_set_removed_nec.size());

  CandidateSetContainer gar_pattern_updated_candidate_set_removed_nec;
  for (const auto& pattern_updated_candidate_removed_nec
                 : pattern_updated_candidate_set_removed_nec) {
    auto pattern_vertex_handle = pattern_updated_candidate_removed_nec.first;
    assert(pattern_vertex_handle);
    auto gar_pattern_vertex_handle
       = gar_pattern.FindVertex(pattern_vertex_handle->id());
    assert(gar_pattern_vertex_handle);
    auto [gar_pattern_updated_candidate_set_removed_nec_it,
          gar_pattern_updated_candidate_set_removed_nec_ret]
         = gar_pattern_updated_candidate_set_removed_nec.emplace(gar_pattern_vertex_handle,
                                            pattern_updated_candidate_removed_nec.second);
    assert(gar_pattern_updated_candidate_set_removed_nec_ret);
  }
  assert(gar_pattern_updated_candidate_set_removed_nec.size()
          == pattern_updated_candidate_set_removed_nec.size());


  int64_t x_supp_increase_counter = 0, x_supp_decrease_counter = 0;
  bool no_increase = false, no_decrease = false;
  if (gar_pattern_added_edge) {
    // further refine the candidate set for more efficient processing
    
    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                         origin_data_graph, 
                                            gar_pattern_origin_candidate_set)) {
      no_decrease = true;
      //std::cout << "no decrease" << std::endl;
    } else {
      gar_pattern_origin_candidate_set_removed_nec = gar_pattern_origin_candidate_set;
      GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                      DataGraphType>(
                              gar_pattern_origin_candidate_set_removed_nec,
                              origin_data_graph_nec);
    }

    if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                         updated_data_graph, 
                                          gar_pattern_updated_candidate_set)) {
      no_increase = true;
      //std::cout << "no increase " << std::endl;
    } else {
      gar_pattern_updated_candidate_set_removed_nec = gar_pattern_updated_candidate_set;
      GUNDAM::RemoveDuplicateCandidate<GraphPatternType, 
                                      DataGraphType>(
                              gar_pattern_updated_candidate_set_removed_nec,
                              updated_data_graph_nec);
    }
  }


  uint64_t x_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter,
           &supp_bound,
           &kConsiderConfidenceBound](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || x_supp_counter < supp_bound);
    x_supp_counter++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (x_supp_counter == supp_bound) {
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  uint64_t xy_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, 
                                            DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter,
           &supp_bound,
           &kConsiderConfidenceBound](
                     const GUNDAM::Match<GraphPatternType,
                                            DataGraphType>& match) {
    assert(kConsiderConfidenceBound
        || xy_supp_counter < supp_bound);
    xy_supp_counter++;
    if (kConsiderConfidenceBound) {
      // need to satisfy the confidence bound, cannot
      // stop when reach supp_bound
      return true;
    }
    // can stop when reached the support bound
    if (xy_supp_counter == supp_bound) {
      // has reached supp_bound
      // stop matching
      return false;
    }
    return true;
  };

  const auto kRhsLiteralInfo
         = (*gar.y_literal_set()
                .begin())->info();

  auto rhs_literal_x_id = kRhsLiteralInfo.x_id();
  auto rhs_literal_y_id = kRhsLiteralInfo.y_id();

  //std::cout << "start to check deleted edges" << std::endl;

  auto begin_time = std::time(NULL);

  std::unordered_map<DataGraphVertexIDType,
                      std::unordered_set<DataGraphVertexIDType>>
      deletion_affected_pivots, all_pivots;


  if (pivot_file_name.length() != 0) {
    std::ifstream pivot_file(pivot_file_name);
    if (!pivot_file.is_open()) {
      std::cout << "can not open input file " << pivot_file_name << std::endl;
      pivot_file_name.clear();
    } else {
      uint64_t data_x_pivot, data_y_pivot, affected_status;
      while (pivot_file >> data_x_pivot >> data_y_pivot >> affected_status) {
        all_pivots[data_x_pivot].insert(data_y_pivot);
        if (affected_status == 1) {
          deletion_affected_pivots[data_x_pivot].insert(data_y_pivot);
        }
      }
      pivot_file.close();
    }
  }
 
  std::map<GraphPatternVertexHandleType,
          std::set<DataGraphVertexHandleType>> origin_affected_area, updated_affected_area;
  std::map<DataGraphVertexHandleType, std::set<DataGraphVertexHandleType>>
        decrease_support_candidate_set;
  //if (pivot_file_name.length() == 0) {
    // get affected pivoted nodes on origin_graph
  for (auto pattern_vertex_it = gar_pattern.VertexBegin();
          !pattern_vertex_it.IsDone();
            pattern_vertex_it++) {

    for (auto pattern_out_edge_it = pattern_vertex_it->OutEdgeBegin();
            !pattern_out_edge_it.IsDone();
              pattern_out_edge_it++) {

      auto pattern_src_handle = pattern_out_edge_it->src_handle();
      auto pattern_dst_handle = pattern_out_edge_it->dst_handle();
      auto pattern_edge_label = pattern_out_edge_it->label();
      auto pattern_dst_label = pattern_dst_handle->label();

      //std::cout << "candidate set size " << gar_pattern_origin_candidate_set[pattern_src_handle].size() << std::endl;

      for (auto data_src_vertex_handle : gar_pattern_origin_candidate_set[pattern_src_handle]) {
        bool has_at_least_one_candidate = false;
        std::vector<DataGraphVertexHandleType> inc_dst_vertex_candidate;

        const auto &pattern_x_handle = gar.pattern().FindVertex(rhs_literal_x_id);
        const auto &pattern_y_handle = gar.pattern().FindVertex(rhs_literal_y_id);
        for (auto data_out_edge_it = data_src_vertex_handle->OutEdgeBegin(pattern_edge_label);
                !data_out_edge_it.IsDone();
                  data_out_edge_it++) {
          auto data_dst_handle = data_out_edge_it->dst_handle();
          //std::cout << "to check an edge " << std::endl;
          const auto &attribute_handle = data_out_edge_it->FindAttribute("inc_state");
          if ((data_dst_handle->label() == pattern_dst_label)
                && data_out_edge_it->template const_attribute<int>("inc_state") == 0) {
            //std::cout << "an edge deleted from " << data_src_vertex_handle->id() << " to " << data_dst_handle->id() << std::endl;
            has_at_least_one_candidate = true;
            inc_dst_vertex_candidate.emplace_back(data_dst_handle);
          }
        }
        if (!has_at_least_one_candidate) {
          continue;
        }
        std::vector<DataGraphVertexHandleType> inc_src_vertex_candidate;
        inc_src_vertex_candidate.emplace_back(data_src_vertex_handle);

        auto current_gar_pattern_origin_candidate_set = gar_pattern_origin_candidate_set;
        // swap with candidate for entire graph
        current_gar_pattern_origin_candidate_set[pattern_src_handle].swap(inc_src_vertex_candidate);
        current_gar_pattern_origin_candidate_set[pattern_dst_handle].swap(inc_dst_vertex_candidate);


        if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                          origin_data_graph, 
                                            current_gar_pattern_origin_candidate_set)) {

          //origin_affected_area;
          for (auto &[pattern_node_handle, origin_graph_node_handle_set]
                      : current_gar_pattern_origin_candidate_set) {
            if ((pattern_node_handle == pattern_x_handle) || (pattern_node_handle == pattern_y_handle)) {
              for (auto &origin_graph_node_handle : origin_graph_node_handle_set) {
                origin_affected_area[pattern_node_handle].insert(origin_graph_node_handle);
              }
            }
          }
        }
      }
    }
  }


  std::map<DataGraphVertexHandleType, std::set<DataGraphVertexHandleType>>
      increase_support_candidate_set;


// get affected pivoted nodes on updated graph
  begin_time = std::time(NULL);
  for (auto pattern_vertex_it = gar_pattern.VertexBegin();
           !pattern_vertex_it.IsDone();
            pattern_vertex_it++) {

    const auto &pattern_x_handle = gar.pattern().FindVertex(rhs_literal_x_id);
    const auto &pattern_y_handle = gar.pattern().FindVertex(rhs_literal_y_id);

    for (auto pattern_out_edge_it = pattern_vertex_it->OutEdgeBegin();
             !pattern_out_edge_it.IsDone();
              pattern_out_edge_it++) {

      auto pattern_src_handle = pattern_out_edge_it->src_handle();
      auto pattern_dst_handle = pattern_out_edge_it->dst_handle();
      auto pattern_edge_label = pattern_out_edge_it->label();
      auto pattern_dst_label = pattern_dst_handle->label();

      for (auto data_src_vertex_handle : gar_pattern_updated_candidate_set[pattern_src_handle]) {
        bool has_at_least_one_candidate = false;
        std::vector<DataGraphVertexHandleType> inc_dst_vertex_candidate;
        for (auto data_out_edge_it = data_src_vertex_handle->OutEdgeBegin(pattern_edge_label);
                 !data_out_edge_it.IsDone();
                  data_out_edge_it++) {
          auto data_dst_handle = data_out_edge_it->dst_handle();
          if ((data_dst_handle->label() == pattern_dst_label)
                && data_out_edge_it->template const_attribute<int>("inc_state") == 2) {
            has_at_least_one_candidate = true;
            inc_dst_vertex_candidate.emplace_back(data_dst_handle);
          }
        }
        if (!has_at_least_one_candidate) {
          continue;
        }

        std::vector<DataGraphVertexHandleType> inc_src_vertex_candidate;
        inc_src_vertex_candidate.emplace_back(data_src_vertex_handle);

        auto current_gar_pattern_updated_candidate_set = gar_pattern_updated_candidate_set;
        current_gar_pattern_updated_candidate_set[pattern_src_handle].swap(inc_src_vertex_candidate);
        current_gar_pattern_updated_candidate_set[pattern_dst_handle].swap(inc_dst_vertex_candidate);


        if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                        updated_data_graph, 
                                           current_gar_pattern_updated_candidate_set)) {
          for (auto &[pattern_node_handle, updated_graph_node_handle_set]
                      : current_gar_pattern_updated_candidate_set) {
            if ((pattern_node_handle == pattern_x_handle) || (pattern_node_handle == pattern_y_handle)) {
              for (auto &updated_graph_node_handle : updated_graph_node_handle_set) {
                updated_affected_area[pattern_node_handle].insert(updated_graph_node_handle);
              }
            }
          }
        }

      }
    }
  }

  for (auto &[pattern_node_handle, origin_node_handle_set] : origin_affected_area) {
    for (auto &origin_node_handle : origin_node_handle_set) {
      auto data_graph_node_id = origin_node_handle->id();
      if (id_to_origin_updated_vertex_handle.find(data_graph_node_id)
              == id_to_origin_updated_vertex_handle.end()) {
        continue;
      }
      auto updated_graph_node_handle
              = id_to_origin_updated_vertex_handle[data_graph_node_id].second;
      if (!updated_graph_node_handle) {
        continue;
      }
      updated_affected_area[pattern_node_handle].insert(updated_graph_node_handle);
    }
  }

  for (auto &[pattern_node_handle, updated_node_handle_set] : updated_affected_area) {
    for (auto &updated_node_handle : updated_node_handle_set) {
      auto data_graph_node_id = updated_node_handle->id();
      if (id_to_origin_updated_vertex_handle.find(data_graph_node_id)
              == id_to_origin_updated_vertex_handle.end()) {
        continue;
      }
      auto origin_graph_node_handle
              = id_to_origin_updated_vertex_handle[data_graph_node_id].first;
      if (!origin_graph_node_handle) {
        continue;
      }
      origin_affected_area[pattern_node_handle].insert(origin_graph_node_handle);
    }
  }

  for (auto &[pattern_node_handle, origin_node_handle_set] : origin_affected_area) {
    std::vector<DataGraphVertexHandleType>
        origin_node_handle_vec(origin_node_handle_set.begin(), origin_node_handle_set.end());
    gar_pattern_origin_candidate_set[pattern_node_handle].swap(origin_node_handle_vec);
  }

  for (auto &[pattern_node_handle, updated_node_handle_set] : updated_affected_area) {
    std::vector<DataGraphVertexHandleType>
        updated_node_handle_vec(updated_node_handle_set.begin(), updated_node_handle_set.end());
    gar_pattern_updated_candidate_set[pattern_node_handle].swap(updated_node_handle_vec);
  }

  int64_t origin_affected_area_supp = 0, updated_affected_area_supp = 0;;

  if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                origin_data_graph, 
                                                gar_pattern_origin_candidate_set)) {
    GUNDAM::Match<GraphPatternType,
                    DataGraphType> partial_match_state;
    x_supp_counter = 0;
    xy_supp_counter = 0;
    auto [temp_x_supp, temp_xy_supp] = 
                gar::RhsGarSupp(gar, origin_data_graph, partial_match_state,
                gar_pattern_origin_candidate_set,
                satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                time_limit_per_supp,
                time_limit_per_supp);
    origin_affected_area_supp = temp_xy_supp;
    //std::cout << "origin affected area supp " << origin_affected_area_supp << std::endl;
  }

  if (GUNDAM::_dp_iso_using_match::RefineCandidateSet(gar_pattern, 
                                                updated_data_graph, 
                                                gar_pattern_updated_candidate_set)) {
    GUNDAM::Match<GraphPatternType,
                    DataGraphType> partial_match_state;
    x_supp_counter = 0;
    xy_supp_counter = 0;
    auto [temp_x_supp, temp_xy_supp] = 
                gar::RhsGarSupp(gar, updated_data_graph, partial_match_state,
                gar_pattern_updated_candidate_set,
                satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                time_limit_per_supp,
                time_limit_per_supp);
    updated_affected_area_supp = temp_xy_supp;
  }

  estimated_support = support + (int64_t)(probability * (double)((int64_t)updated_affected_area_supp - (int64_t)origin_affected_area_supp));

  if (support >= supp_bound) {
    int64_t origin_upper_bound = gar::EstimateUpperBound(probability, origin_affected_area_supp, kThreshold);
    int64_t updated_lower_bound = gar::EstimateLowerBound(probability, updated_affected_area_supp, kThreshold);

    if (support - origin_upper_bound + updated_lower_bound >= supp_bound) {
      return true;
    }
    int64_t origin_lower_bound = gar::EstimateLowerBound(probability, origin_affected_area_supp, kThreshold);
    int64_t updated_upper_bound = gar::EstimateUpperBound(probability, updated_affected_area_supp, kThreshold);
    if (support - origin_lower_bound + updated_upper_bound < supp_bound) {
      return true;
    }
    return false;
  } else {
    int64_t origin_lower_bound = gar::EstimateLowerBound(probability, origin_affected_area_supp, kThreshold);
    int64_t updated_upper_bound = gar::EstimateUpperBound(probability, updated_affected_area_supp, kThreshold);

    if (support - origin_lower_bound + updated_upper_bound < supp_bound) {
      return true;
    }
    int64_t origin_upper_bound = gar::EstimateUpperBound(probability, origin_affected_area_supp, kThreshold);
    int64_t updated_lower_bound = gar::EstimateLowerBound(probability, updated_affected_area_supp, kThreshold);
    if (support - origin_upper_bound + updated_lower_bound >= supp_bound) {
      return true;
    }
    return false;
  }
  return true;
}













}  // namespace _gar_discover

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_HAS_MATCH_
