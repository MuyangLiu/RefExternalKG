#ifndef GAR_SUPP_H_
#define GAR_SUPP_H_

#include <map>
#include <unistd.h>
#include <cmath>

#include "gar/gar.h"
#include "gar/same_gar.h"

#include "gundam/algorithm/match_using_match.h"

#include "gundam/match/match.h"

#include "gundam/type_getter/vertex_handle.h"
#include "gundam/type_getter/edge_handle.h"
#include "gundam/tool/operator/preserve.h"
#include "gundam/tool/sub_graph_of.h"

#include "gundam/io/csvgraph.h"

namespace gar{
namespace _gar_supp{

template <typename LiteralContainerType,
          typename          VertexIDType>
void AddVertexIDTo(LiteralContainerType& literal_container,
                  std::set<VertexIDType>& vertex_id_set){
  for (const auto& literal_ptr : literal_container) {
    auto literal_info = literal_ptr->info();
    switch (literal_info.literal_type()){
    case gar::LiteralType::kAttrValueLiteral:
      // x.A
    case gar::LiteralType::kConstantLiteral:
      // x.A == c
      vertex_id_set.emplace(literal_info.x_id());
      continue;

    case gar::LiteralType::kEdgeLiteral:
      // x_id_ - edge_label_ -> y_id_
    case gar::LiteralType::kVariableLiteral:
      // x.A == y.B
    case gar::LiteralType::kMlLiteral:
      // x_id_ - (module_url_: edge_label_) -> y_id_
      vertex_id_set.emplace(literal_info.x_id());
      vertex_id_set.emplace(literal_info.y_id());
      continue;
    case gar::LiteralType::kKeyLiteral:
      // x_id_ == y_id_
      vertex_id_set.emplace(literal_info.x_id());
      vertex_id_set.emplace(literal_info.y_id());
      continue;
    default:
      // unknown literal type
      std::cout << "unknown literal type" << std::endl;
      assert(false);
      break;
    }
  }
  return;
}

/* ############################################################# *
 * ##  this function is utilized to grep the sub-pattern      ## *
 * ##  that are contained in the literals and user-specified  ## *
 * ##  pivot set                                              ## *     
 * ##  for example, consider the following example gar:       ## *
 * ##     pattern:  0 -> 1 -> 2 -> 3                          ## *
 * ##           X:  0.A = a0                                  ## *
 * ##           Y:  2.A = a1                                  ## *
 * ##   when pivoted_vertex_set is empty                      ## *
 * ##     rhs_only == true                                    ## *
 * ##       return 2  (isolated vertex)                       ## *
 * ##     rhs_only == false                                   ## *
 * ##       return 0   2  (two isolated vertexex)             ## *
 * ##                                                         ## *
 * ##   when pivoted_vertex_set == {1}                        ## *
 * ##     rhs_only == true                                    ## *
 * ##       return 1 -> 2                                     ## *
 * ##     rhs_only == false                                   ## *
 * ##       return 0 -> 1 -> 2                                ## *
 * ############################################################# */
template <typename GraphPatternType,
          typename    DataGraphType>
inline GraphPatternType LiteralPattern(
  GraphAssociationRule<GraphPatternType,
                          DataGraphType>& gar, bool rhs_only,
  const std::vector<typename 
  GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set
      = std::vector<typename 
  GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using VertexIDType = typename GUNDAM::VertexID<GraphPatternType>::type;

  std::set<VertexIDType> preserve_vertex_id_set;

  AddVertexIDTo(gar.y_literal_set(), preserve_vertex_id_set);
  if (!rhs_only) {
    AddVertexIDTo(gar.x_literal_set(), preserve_vertex_id_set);
  }
  
  for (const auto& vertex_handle : pivoted_vertex_set) {
    preserve_vertex_id_set.emplace(vertex_handle->id());
  }

  std::set<VertexIDType> erase_vertex_id_set;
  for (auto vertex_it = gar.pattern().VertexBegin();
           !vertex_it.IsDone();
            vertex_it++) {
    if (preserve_vertex_id_set.find(vertex_it->id())
     != preserve_vertex_id_set.end()){
      // contained in preserve_vertex_id_set
      continue;
    }
    erase_vertex_id_set.emplace(vertex_it->id());
  }
  auto literal_pattern(gar.pattern());
  for (const auto& erase_vertex_id : erase_vertex_id_set){
    literal_pattern.EraseVertex(erase_vertex_id);
  }
  return literal_pattern;
}

template <typename GraphPatternType,
          typename    DataGraphType>
inline GraphAssociationRule<GraphPatternType,
                               DataGraphType> 
RhsGar(GraphAssociationRule<GraphPatternType,
                               DataGraphType>& gar,
       const std::vector<typename 
       GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set
           = std::vector<typename 
       GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
    pivoted_rhs_gar(_gar_supp::LiteralPattern(gar, true, pivoted_vertex_set));
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_rhs_pattern(gar.pattern(), 
                   pivoted_rhs_gar.pattern(), "same_id_map");
  assert(gar_pattern_to_rhs_pattern.size()
             == pivoted_rhs_gar.pattern().CountVertex());
  for (const auto& y_literal_ptr : gar.y_literal_set()) {
    #ifndef NDEBUG
    std::vector<GraphPatternVertexHandle> y_literal_vertex_handle_set;
    y_literal_ptr->CalPivot(y_literal_vertex_handle_set);
    assert(y_literal_vertex_handle_set.size() == 1
        || y_literal_vertex_handle_set.size() == 2);
    for (const auto& y_handle : y_literal_vertex_handle_set) {
      assert(gar_pattern_to_rhs_pattern.HasMap(y_handle));
    }
    #endif // NDEBUG
    pivoted_rhs_gar.AddY(y_literal_ptr->info());
  }
  for (const auto& x_literal_ptr : gar.x_literal_set()) {
    std::vector<GraphPatternVertexHandle> x_literal_vertex_handle_set;
    x_literal_ptr->CalPivot(x_literal_vertex_handle_set);
    assert(x_literal_vertex_handle_set.size() == 1
        || x_literal_vertex_handle_set.size() == 2);
    bool x_literal_mapped = true;
    for (const auto& x_handle : x_literal_vertex_handle_set) {
      if (!gar_pattern_to_rhs_pattern.HasMap(x_handle)) {
        x_literal_mapped = false;
        break;
      }
    }
    if (!x_literal_mapped) {
      continue;
    }
    // add all x literals that are contained 
    pivoted_rhs_gar.AddX(x_literal_ptr->info());
  }
  return pivoted_rhs_gar;
}

template <typename GraphPatternType,
          typename    DataGraphType>
inline GraphAssociationRule<GraphPatternType,
                               DataGraphType> 
LiteralGar(GraphAssociationRule<GraphPatternType,
                                   DataGraphType>& gar,
           const std::vector<typename 
           GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
               = std::vector<typename 
           GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
    literal_gar(_gar_supp::LiteralPattern(gar, false, pivoted_vertex_set));
  GUNDAM::Match<GraphPatternType, 
                GraphPatternType> 
    gar_pattern_to_literal_pattern(gar.pattern(), 
                           literal_gar.pattern(),
                           "same_id_map");
  assert(gar_pattern_to_literal_pattern.size()
        == literal_gar.pattern().CountVertex());
  for (const auto& y_literal_ptr   
             : gar.y_literal_set()){
    #ifndef NDEBUG
    std::vector<GraphPatternVertexHandle> y_literal_vertex_handle_set;
    y_literal_ptr->CalPivot(y_literal_vertex_handle_set);
    assert(y_literal_vertex_handle_set.size() == 1
        || y_literal_vertex_handle_set.size() == 2);
    for (const auto& y_handle : y_literal_vertex_handle_set) {
      assert(gar_pattern_to_literal_pattern.HasMap(y_handle));
    }
    #endif // NDEBUG
    literal_gar.AddY(y_literal_ptr->info());
  }
  for (const auto& x_literal_ptr 
             : gar.x_literal_set()){
    std::vector<GraphPatternVertexHandle> x_literal_vertex_handle_set;
    x_literal_ptr->CalPivot(x_literal_vertex_handle_set);
    assert(x_literal_vertex_handle_set.size() == 1
        || x_literal_vertex_handle_set.size() == 2);
    bool x_literal_mapped = true;
    for (const auto& x_handle : x_literal_vertex_handle_set) {
      if (!gar_pattern_to_literal_pattern.HasMap(x_handle)){
        x_literal_mapped = false;
        break;
      }
    }
    if (!x_literal_mapped){
      continue;
    }
    literal_gar.AddX(x_literal_ptr->info());
  }
  return std::move(literal_gar);
}

} // namespace _gar_supp

/*  ################################################  *
 *  ##                                            ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##          pivoted_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##      pivoted_rhs_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##  pivoted_literal_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##              gar_pattern_match_callback    ##  *
 *  ##                                            ##  *
 *  ################################################  */

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                                     const bool store_match,
     std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>>& match_vec,
     double& probability,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  assert(time_limit == -1.0
      || time_limit > 0);

  assert(time_limit_per_supp == -1.0
      || time_limit_per_supp > 0);

  assert(pivoted_vertex_set.empty() 
      || pivoted_vertex_set.size() > gar_pattern_to_data_graph_partial_match.size());

  if (// time_limit is set but time_limit_per_supp is not set
      (time_limit_per_supp == -1.0 && time_limit != -1.0)
      // both time_limit and time_limit_per_supp are set
      // but there are time_limit_per_supp > time_limit
   || (time_limit_per_supp != -1.0 && time_limit != -1.0
    && time_limit_per_supp > time_limit)) {
    time_limit_per_supp = time_limit;
  }

  auto& gar_pattern = gar.pattern();

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;

  /* ############################################################## *
   * ##   obtain pattern contained in all literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
        pivoted_literal_gar = _gar_supp::LiteralGar(gar, pivoted_vertex_set);
  auto& pivoted_literal_pattern = pivoted_literal_gar.pattern();

  /* ############################################################## *
   * ##   obtain pattern contained in rhs literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
         pivoted_rhs_gar = _gar_supp::RhsGar(gar, pivoted_vertex_set);
  auto&  pivoted_rhs_pattern = pivoted_rhs_gar.pattern();
  assert(pivoted_literal_pattern.CountVertex()
          >= pivoted_rhs_pattern.CountVertex());
  assert(GUNDAM::SubGraphOf(pivoted_rhs_pattern,
                        pivoted_literal_pattern));

  /* ############################################# *
   * ##   obtain pattern contained in pivotes   ## *
   * ############################################# */
  GraphPatternType pivoted_pattern = pivoted_vertex_set.empty()?
                                     pivoted_rhs_pattern
                                   : GUNDAM::PreserveVertexSet(gar.pattern(),
                                     pivoted_vertex_set);
  assert(pivoted_rhs_pattern.CountVertex() 
          >= pivoted_pattern.CountVertex());
  assert(pivoted_vertex_set.empty()
      || pivoted_pattern.CountVertex() 
      == pivoted_vertex_set.size());
  assert(GUNDAM::SubGraphOf(pivoted_pattern,
                        pivoted_rhs_pattern));
  #ifndef NDEBUG
  for (auto map_it = gar_pattern_to_data_graph_partial_match.MapBegin();
           !map_it.IsDone();
            map_it++) {
    assert(pivoted_pattern.FindVertex(map_it->src_handle()->id()));
  }
  #endif // NDEBUG

  /* ########################################## *
   * ##  parameters to be used in callbacks  ## *
   * ########################################## */
  const bool pivoted_pattern_same_as_pivoted_rhs_pattern 
               = GUNDAM::SamePattern(pivoted_pattern,
                                     pivoted_rhs_pattern);

  const bool pivoted_rhs_gar_same_as_pivoted_literal_gar 
                           = SameGar(pivoted_rhs_gar,
                                     pivoted_literal_gar);

  const bool pivoted_literal_pattern_same_as_gar_pattern
           = GUNDAM::SamePattern(pivoted_literal_pattern,
                                             gar_pattern);

  /* ############################ *
   * ##  from pivoted_pattern  ## *
   * ############################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match(
    pivoted_rhs_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_pattern_partial_match(
    gar_pattern,   pivoted_pattern, "same_id_map");

  /* ################################ *
   * ##  from pivoted_rhs_pattern  ## *
   * ################################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_rhs_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_rhs_pattern_partial_match(
    gar_pattern,   pivoted_rhs_pattern, "same_id_map");

  /* #################################### *
   * ##  from pivoted_literal_pattern  ## *
   * #################################### */
  GUNDAM::Match<GraphPatternType, GraphPatternType>
    gar_pattern_to_pivoted_literal_pattern_partial_match(
    gar_pattern,   pivoted_literal_pattern, "same_id_map");

  CandidateSetContainer pivoted_pattern_candidate_set,
                    pivoted_rhs_pattern_candidate_set,
                pivoted_literal_pattern_candidate_set;

  // process pivoted_literal_pattern_candidate_set
  //         and pivoted_rhs_pattern_candidate_set
  for (const auto& [gar_pattern_candidate_ptr,
                    gar_pattern_candidate]
                  : gar_pattern_candidate_set_for_vertexes_contained_in_literals) {
    // should not be null
    assert(gar_pattern_candidate_ptr);
    const auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
    if (kPivotedLiteralPatternPtr) {
      // this vertex is contained in pivoted_literal_pattern
      auto [ pivoted_literal_pattern_candidate_set_it,
             pivoted_literal_pattern_candidate_set_ret ]
           = pivoted_literal_pattern_candidate_set.emplace(
            kPivotedLiteralPatternPtr,
            gar_pattern_candidate);
      // should added successfully
      assert(pivoted_literal_pattern_candidate_set_ret);
      const auto kPivotedRhsPatternPtr
                = pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id());
      if (kPivotedRhsPatternPtr){
        // this vertex is contained in pivoted_rhs_pattern
        auto [ pivoted_rhs_pattern_candidate_set_it,
               pivoted_rhs_pattern_candidate_set_ret ]
             = pivoted_rhs_pattern_candidate_set.emplace(
              kPivotedRhsPatternPtr,
                gar_pattern_candidate);
        // should added successfully
        assert(pivoted_rhs_pattern_candidate_set_ret);
        const auto kPivotedPatternPtr
                  = pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedPatternPtr){
          // this vertex is contained in pivoted_rhs_pattern
          auto [ pivoted_pattern_candidate_set_it,
                 pivoted_pattern_candidate_set_ret ]
               = pivoted_pattern_candidate_set.emplace(
                 kPivotedPatternPtr,
                  gar_pattern_candidate);
          // should added successfully
          assert(pivoted_pattern_candidate_set_ret);
        }
      }
      #ifndef NDEBUG
      else {
        // should not be contained in pivoted_rhs_pattern
        assert(!pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      }
      #endif // NDEBUG
    }
    #ifndef NDEBUG
    else {
      // should not be contained in pivoted_rhs_pattern and pivoted_pattern
      assert(!pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      assert(    !pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
    }
    #endif // NDEBUG
  }

  bool has_satisfy_x       = false,
       has_satisfy_x_not_y = false,
       has_satisfy_x_and_y = false;

  typename DataGraphType::VertexCounterType x_supp = 0,
                                           xy_supp = 0;

  auto pivoted_rhs_pattern_x_id = (pivoted_rhs_pattern.VertexBegin())->id();
  auto pivoted_rhs_pattern_y_id = (++pivoted_rhs_pattern.VertexBegin())->id();

  uint64_t rhs_gar_supp = 0;

  std::map<typename GUNDAM::VertexID<DataGraphType>::type,
          std::set<typename GUNDAM::VertexID<DataGraphType>::type>> rhs_gar_supp_map;

  /* ################################################ *
   * ##  callbacks from gar_pattern to data graph  ## *
   * ################################################ */
  std::function<bool(const MatchType&)> 
    gar_pattern_prune_callback
    = [](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // prune nothing, continue the matching
    return false;
  };
  // based on the partial match from literal_kPattern to data graph
  // complete the match of the entire gar_pattern 
  std::function<bool(const MatchType&)> 
    gar_pattern_match_callback 
    = [&has_satisfy_x,
       &store_match,
       &match_vec](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // only need to complete the match, find one match that 
    // satisfy all literal is enough
    has_satisfy_x = true;

    if (store_match) {
      match_vec.emplace_back(gar_pattern_to_data_graph_match);
    }
    // terminate matching
    return false;
  };

  /* ############################################################ *
   * ##  callbacks from pivoted_literal_pattern to data graph  ## *
   * ############################################################ */
  // if the current pivoted_literal_pattern_to_data_graph_match
  // has already violate a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : pivoted_literal_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_match)){
        // this literal is not considered in the current
        // pivoted_literal_pattern_to_data_graph_match,
        // this match cannot be pruned
        //
        // move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    if (!pivoted_literal_pattern_same_as_gar_pattern) {
      // pivoted_rhs_gar is not the same as pivoted_literal_gar
      // first match the pivoted_rhs_literal, then match
      // the pivoted_literal_gar
      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
        = pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_pivoted_literal_pattern_partial_match);
      assert(!has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
      if (has_satisfy_x) {
        // has found satisfy X, terminate matching
        return false;
      }
      // continue matching
      return true;
    }
    assert(!has_satisfy_x);
    // gar_pattern of literal_gar is the same as the
    // entire gar_pattern, does not need further match
    has_satisfy_x = true;

    if (store_match) {
      GUNDAM::Match<GraphPatternType, DataGraphType>
                      gar_pattern_to_data_graph_entire_match 
        = pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_pivoted_literal_pattern_partial_match);
      match_vec.emplace_back(gar_pattern_to_data_graph_entire_match);
    }
    // has found satisfy X, terminate matching
    return false;
  };

  /* ######################################################## *
   * ##  callbacks from pivoted_rhs_pattern to data graph  ## *
   * ######################################################## */
  // if the current match_state has already violate
  // a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_prune_callback
    = [&](const MatchType&
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : pivoted_rhs_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match)) {
        // this literal is not considered in the current
        // pivoted_rhs_pattern_to_data_graph_match, 
        // this match cannot be pruned move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_match_callback
    = [&](const MatchType& 
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    assert(!has_satisfy_x);
    assert(!has_satisfy_x_not_y);
    /* ############################################## *
     * ##  first check whether satisfy y literals  ## *
     * ############################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }

    if (!pivoted_rhs_gar_same_as_pivoted_literal_gar) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar

      /* ####################################################### *
       * ##  get the paritial match from pivoted_rhs_pattern  ## *
       * ##  to data_graph                                    ## *
       * ####################################################### */
      GUNDAM::Match<GraphPatternType, DataGraphType>
      pivoted_literal_pattern_to_data_graph_partial_match
        = pivoted_rhs_pattern_to_data_graph_match(
      pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match);
      
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph, 
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set,
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    }
    else {
      // the pivoted_rhs_gar is the same as pivoted_literal_gar
      // directly match literal_gar to data_graph
      if (!pivoted_literal_pattern_same_as_gar_pattern) {
        /* ####################################################### *
         * ##  get the paritial match from pivoted_rhs_pattern  ## *
         * ##  to data_graph                                    ## *
         * ####################################################### */
        GUNDAM::Match<GraphPatternType, DataGraphType>
                  gar_pattern_to_data_graph_partial_match
        = pivoted_rhs_pattern_to_data_graph_match(
                  gar_pattern_to_pivoted_rhs_pattern_partial_match);

        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                gar_pattern,   data_graph, 
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
      }
      else {
        if (store_match) {
          GUNDAM::Match<GraphPatternType, DataGraphType>
                          gar_pattern_to_data_graph_entire_match 
            = pivoted_rhs_pattern_to_data_graph_match(
                          gar_pattern_to_pivoted_rhs_pattern_partial_match);
          match_vec.emplace_back(gar_pattern_to_data_graph_entire_match);
        }

        // pivoted_rhs_gar_pattern is the same as gar_pattern
        has_satisfy_x = true;
      }
    }
    assert(!has_satisfy_x_not_y);
    assert(!has_satisfy_x_and_y);
    if (!has_satisfy_x) {
      // continue matching
      return true;
    }
    // satisfy x
    // reset parameter
    has_satisfy_x = false;
    if (satisify_y_literals) {
      // satisfy x and y
      assert(!has_satisfy_x_not_y);
      has_satisfy_x_and_y = true;
      // continue matching
      return true;
    }
    // satisfy x but not y
    has_satisfy_x_and_y = false;
    has_satisfy_x_not_y =  true;
    // no longer needs matching
    return false;
  };

  /* #################################################### *
   * ##  callbacks from pivoted_pattern to data graph  ## *
   * #################################################### */
  std::function<bool(const MatchType&)>
    pivoted_pattern_prune_callback
    = [](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    // prune nothing
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_pattern_match_callback
    = [&](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    /* ######################################## *
     * ##  get a match for pivot vertex set  ## *
     * ######################################## */
    has_satisfy_x       = false;
    has_satisfy_x_not_y = false;
    has_satisfy_x_and_y = false;
    
    /* ####################################################### *
     * ##  get the paritial match from pivoted_rhs_pattern  ## *
     * ##  to data_graph                                    ## *
     * ####################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_rhs_pattern_to_data_graph_partial_match
      = pivoted_pattern_to_data_graph_match(
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match);

    if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
      std::cout << "should not reach here" << std::endl;
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_rhs_pattern,   data_graph, 
                              pivoted_rhs_pattern_to_data_graph_partial_match,
                              pivoted_rhs_pattern_candidate_set,
                              pivoted_rhs_pattern_prune_callback,
                              pivoted_rhs_pattern_match_callback,
                              time_limit_per_supp);
      assert(!has_satisfy_x);
      assert(!has_satisfy_x_and_y 
          || !has_satisfy_x_not_y); // cannot be both true
      has_satisfy_x = has_satisfy_x_and_y 
                   || has_satisfy_x_not_y;
      if (!has_satisfy_x) {
        return true;
      }
      x_supp++;
      const bool kSatisfyXSuppCallbackRet
                = satisfy_x_supp_callback(pivoted_pattern_to_data_graph_match);
      if (has_satisfy_x_not_y) {
        return kSatisfyXSuppCallbackRet;
      }
      assert(has_satisfy_x_and_y);
      xy_supp++;
      const bool kSatisfyXYSuppCallbackRet
                = satisfy_xy_supp_callback(pivoted_pattern_to_data_graph_match);
      return kSatisfyXSuppCallbackRet
         && kSatisfyXYSuppCallbackRet;
    }
    /* ########################################################### *
     * ##  pivoted_pattern is the same as pivoted_rhs_pattern,  ## *
     * ##  directly whether the match satisfy y literals here   ## *
     * ########################################################### */
    assert(pivoted_rhs_pattern_to_data_graph_partial_match.size()
            == pivoted_pattern_to_data_graph_match.size());

    /* ########################################################### *
     * ##  get the paritial match from pivoted_literal_pattern  ## *
     * ##  to data_graph                                        ## *
     * ########################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
    pivoted_literal_pattern_to_pivoted_pattern_partial_match);


    /* ######################################################## *
    * ##  chech whether the match satisfy rhs literal here  ## *
    * ######################################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_partial_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_partial_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }
    if (satisify_y_literals) {
      auto data_graph_x_handle = pivoted_rhs_pattern_to_data_graph_partial_match.MapTo(pivoted_rhs_pattern_x_id);
      auto data_graph_y_handle = pivoted_rhs_pattern_to_data_graph_partial_match.MapTo(pivoted_rhs_pattern_y_id);
      if (data_graph_x_handle && data_graph_y_handle) {
        auto data_graph_x_id = data_graph_x_handle->id();
        auto data_graph_y_id = data_graph_y_handle->id();
        if (rhs_gar_supp_map[data_graph_x_id].find(data_graph_y_id)
                      == rhs_gar_supp_map[data_graph_x_id].end()) {
          rhs_gar_supp++;
          rhs_gar_supp_map[data_graph_x_id].insert(data_graph_y_id);
        }
      }

    }

    if (!pivoted_rhs_gar_same_as_pivoted_literal_gar) {
      assert( GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_pattern));
      assert(!        SameGar    (pivoted_rhs_gar,     pivoted_literal_gar));
      assert(!GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_literal_pattern));
      
      assert(!has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph, 
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set,
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    }
    else {
      // pivoted_rhs_pattern is the same as pivoted_literal_pattern
      // directly match to pivoted_literal_pattern
      /* ############################################### *
       * ##  get the paritial match from gar_pattern  ## *
       * ##  to data_graph                            ## *
       * ############################################### */
      assert(        SameGar    (pivoted_rhs_gar, pivoted_literal_gar));
      assert(GUNDAM::SamePattern(pivoted_pattern, pivoted_literal_pattern));
      /* ######################################################## *
       * ##  chech whether the match satisfy lhs literal here  ## *
       * ######################################################## */
      bool satisify_x_literals = true;
      for (const auto& x_literal_ptr : pivoted_literal_gar.x_literal_set()) {
        assert(x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_partial_match));
        if (x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_partial_match)) {
          continue;
        }
        satisify_x_literals = false;
        break;
      }
      if (!satisify_x_literals) {
        // continue matching
        assert(!has_satisfy_x);
        return true;
      }
      if (!pivoted_literal_pattern_same_as_gar_pattern) {
        GUNDAM::Match<GraphPatternType, DataGraphType>
              gar_pattern_to_data_graph_partial_match
        = pivoted_pattern_to_data_graph_match(
              gar_pattern_to_pivoted_pattern_partial_match);
        assert(!has_satisfy_x);
        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                gar_pattern,   data_graph, 
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
      }
      else {
        assert(GUNDAM::SamePattern(pivoted_pattern, gar_pattern));
        // pivoted_pattern is the same as gar_pattern, no longer need to match
        assert(!has_satisfy_x);

        if (store_match) {
          GUNDAM::Match<GraphPatternType, DataGraphType>
                          gar_pattern_to_data_graph_entire_match 
            = pivoted_pattern_to_data_graph_match(
                          gar_pattern_to_pivoted_pattern_partial_match);
          match_vec.emplace_back(gar_pattern_to_data_graph_entire_match);
        }

        has_satisfy_x = true;
      }
    }
    assert(!has_satisfy_x_and_y);
    assert(!has_satisfy_x_not_y);

    if (!has_satisfy_x) {
      return true;
    }
    x_supp++;

    GUNDAM::Match<GraphPatternType, DataGraphType>
          gar_pattern_to_data_graph_partial_match
    = pivoted_pattern_to_data_graph_match(
          gar_pattern_to_pivoted_pattern_partial_match);

    const bool kSatisfyXSuppCallbackRet
              = satisfy_x_supp_callback(gar_pattern_to_data_graph_partial_match);

    /* ######################################################## *
    * ##  chech whether the match satisfy rhs literal here  ## *
    * ######################################################## */
    // bool satisify_y_literals = true;
    // for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
    //   assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_partial_match));
    //   if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_partial_match)) {
    //     continue;
    //   }
    //   satisify_y_literals = false;
    //   break;
    // }
    if (!satisify_y_literals) {
      return kSatisfyXSuppCallbackRet;
    }
    // satisfy y literals
    xy_supp++;
    // continue matching
    const bool kSatisfyXYSuppCallbackRet
              = satisfy_xy_supp_callback(gar_pattern_to_data_graph_partial_match);
    return kSatisfyXSuppCallbackRet
       && kSatisfyXYSuppCallbackRet;
  };

  MatchType pivoted_pattern_partial_match;
  for (auto pivoted_pattern_vertex_it = pivoted_pattern.VertexBegin();
           !pivoted_pattern_vertex_it.IsDone();
            pivoted_pattern_vertex_it++) {
    auto target_handle = gar_pattern_to_data_graph_partial_match.MapTo(pivoted_pattern_vertex_it->id());
    if (!target_handle) {
      continue;
    }
    pivoted_pattern_partial_match.AddMap(pivoted_pattern_vertex_it,
                                         target_handle);
  }

  // first match pivoted gar_pattern, then match rhs gar_pattern
  GUNDAM::MatchUsingMatch<match_semantics,
                          GUNDAM::MatchAlgorithm::kDagDp,
                          GUNDAM::MergeNecConfig::kNotMerge>(
                          pivoted_pattern, data_graph,
                          pivoted_pattern_partial_match,
                          pivoted_pattern_candidate_set,
                          pivoted_pattern_prune_callback,
                          pivoted_pattern_match_callback,
                          time_limit);
  assert(x_supp >= xy_supp);
  if (xy_supp == 0 || rhs_gar_supp == 0) {
    probability = 0.0;
  } else {
    probability = (double)xy_supp / (double)rhs_gar_supp;
  }
  // if (xy_supp > 200) {
  //   std::cout << "rhs_gar_supp " << rhs_gar_supp << " xy_supp " << xy_supp << " probability " << probability << std::endl;
  // }
  return std::pair(x_supp, xy_supp);
}



// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  bool store_match = false;
  std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>> match_vec;
  double probability = 0.0;
  return GarSupp(gar, data_graph, gar_pattern_to_data_graph_partial_match,
                gar_pattern_candidate_set_for_vertexes_contained_in_literals,
                gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                store_match, match_vec, probability,
                satisfy_x_supp_callback, satisfy_xy_supp_callback,
                time_limit, time_limit_per_supp, pivoted_vertex_set);
}

// return std::pair<x_supp, xy_supp>
// template <enum GUNDAM::MatchSemantics match_semantics 
//              = GUNDAM::MatchSemantics::kIsomorphism,
//           typename GraphPatternType,
//           typename    DataGraphType>
// std::pair<typename    DataGraphType::VertexCounterType,
//           typename    DataGraphType::VertexCounterType> 
//   GarSupp(GraphAssociationRule<GraphPatternType,
//                                   DataGraphType>& gar,
//              DataGraphType& data_graph,
//   const GUNDAM::Match<GraphPatternType,
//                          DataGraphType>& gar_pattern_to_data_graph_partial_match,
//   const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
//      std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
//   const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
//      std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
//      std::function<bool(const GUNDAM::Match<GraphPatternType,
//                                                DataGraphType>&)>  satisfy_x_supp_callback,
//      std::function<bool(const GUNDAM::Match<GraphPatternType,
//                                                DataGraphType>&)> satisfy_xy_supp_callback,
//             double time_limit,
//             double time_limit_per_supp,
//             const bool store_match,
//             std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>>& match_vec) {

//   std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type> pivoted_vertex_set;
//   return GarSupp<match_semantics>
//                (gar, data_graph, 
//                 gar_pattern_to_data_graph_partial_match,
//                 gar_pattern_candidate_set_for_vertexes_contained_in_literals,
//                 gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
//                 satisfy_x_supp_callback,
//                 satisfy_xy_supp_callback,
//                 time_limit,
//                 time_limit_per_supp,
//                 pivoted_vertex_set,
//                 store_match,
//                 match_vec);
// }

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
                                  DataGraphType&  data_graph,
           const GUNDAM::Match<GraphPatternType,
                                  DataGraphType>& pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  return GarSupp<match_semantics>
                (gar, data_graph, 
                 pattern_to_data_graph_partial_match,
                 pattern_candidate_set,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                 satisfy_xy_supp_callback,
                 time_limit,
                 time_limit_per_supp,
                 pivoted_vertex_set);
}

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
                                  DataGraphType&  data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType> match_state;
  return GarSupp<match_semantics>
                (gar, data_graph, 
                 match_state,
                 pattern_candidate_set,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                 satisfy_xy_supp_callback,
                 time_limit,
                 time_limit_per_supp,
                 pivoted_vertex_set);
}

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
                                  DataGraphType& data_graph,
           const GUNDAM::Match<GraphPatternType,
                                  DataGraphType>& pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
            int64_t  x_supp_limit = -1,
            int64_t xy_supp_limit = -1,
            double     time_limit = -1.0,
            double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){

  using MatchType = GUNDAM::Match<GraphPatternType,
                                     DataGraphType>;

  assert(x_supp_limit == -1 
     || xy_supp_limit == -1 
     || (x_supp_limit >= xy_supp_limit));

  int64_t x_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter,
           &x_supp_limit](const MatchType& match) {
    if (x_supp_limit == -1){
      // x_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(x_supp_counter < x_supp_limit);
    x_supp_counter++;
    if (x_supp_counter == x_supp_limit){
      // has reached x_supp_limit
      // stop matching
      return false;
    }
    return true;
  };

  int64_t xy_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter,
           &xy_supp_limit](const MatchType& match) {
    if (xy_supp_limit == -1) {
      // xy_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(xy_supp_counter < xy_supp_limit);
    xy_supp_counter++;
    if (xy_supp_counter == xy_supp_limit) {
      // has reached xy_supp_limit
      // stop matching
      return false;
    }
    return true;
  };
  
  auto [x_supp, xy_supp] = GarSupp<match_semantics>(
                                   gar, data_graph,
                             pattern_to_data_graph_partial_match,
                             pattern_candidate_set,
                             satisfy_x_supp_callback,
                            satisfy_xy_supp_callback,
                                  time_limit,
                                  time_limit_per_supp,
                                  pivoted_vertex_set);

  assert((( x_supp_limit == -1) || ( x_supp ==  x_supp_counter) && ( x_supp <=  x_supp_limit))
      && ((xy_supp_limit == -1) || (xy_supp == xy_supp_counter) && (xy_supp <= xy_supp_limit)));

  return std::pair(x_supp, xy_supp);
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
     std::function<bool(const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type, 
                                       typename GUNDAM::VertexHandle<   DataGraphType>::type>&)> satisfy_x_supp_callback,
     std::function<bool(const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type, 
                                       typename GUNDAM::VertexHandle<   DataGraphType>::type>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){

  const GUNDAM::Match<GraphPatternType,
                         DataGraphType> pattern_to_data_graph_partial_match;

  return GarSupp<match_semantics>(
                 gar, data_graph,
                 pattern_to_data_graph_partial_match,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                              time_limit,
                              time_limit_per_supp,
                           pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& pattern_to_data_graph_partial_match,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){
  
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using CandidateSetType = std::map<GraphPatternVertexHandle, 
                        std::vector<   DataGraphVertexHandle>>;

  CandidateSetType pattern_candidate_set;
  if (!GUNDAM::_dp_iso::InitCandidateSet<match_semantics>(
                gar.pattern(),
                data_graph,
                pattern_candidate_set)) {
    return std::pair(0, 0);
  }
  if (!GUNDAM::_dp_iso::RefineCandidateSet(
                gar.pattern(),
                data_graph,
                pattern_candidate_set)) {
    return std::pair(0, 0);
  }

  return GarSupp<match_semantics>(
                 gar, data_graph,
                 pattern_to_data_graph_partial_match,
                 pattern_candidate_set,
                 satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                              time_limit,
                              time_limit_per_supp,
                           pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()){

  const GUNDAM::Match<GraphPatternType,
                         DataGraphType> pattern_to_data_graph_partial_match;

  return GarSupp<match_semantics>(
                 gar, data_graph,
                 pattern_to_data_graph_partial_match,
                 satisfy_x_supp_callback,
                satisfy_xy_supp_callback,
                              time_limit,
                              time_limit_per_supp,
                           pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& pattern_candidate_set,
            int64_t  x_supp_limit = -1,
            int64_t xy_supp_limit = -1,
            double     time_limit = -1.0,
            double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  GUNDAM::Match<GraphPatternType,
                   DataGraphType> pattern_to_data_graph_partial_match;
  return GarSupp<match_semantics>(
                 gar, data_graph, pattern_to_data_graph_partial_match, 
                                  pattern_candidate_set, 
                                           x_supp_limit, 
                                          xy_supp_limit,
                                             time_limit,
                                             time_limit_per_supp,
                                          pivoted_vertex_set);
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  GarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
            int64_t  x_supp_limit = -1,
            int64_t xy_supp_limit = -1,
            double     time_limit = -1.0,
            double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;

  auto& pattern  = gar.pattern();

  CandidateSetContainer pattern_candidate_set;
  if (!GUNDAM::_dp_iso::InitCandidateSet<match_semantics>(
                pattern ,
                data_graph,
                pattern_candidate_set)) {
    return std::make_pair(0, 0);
  }
  if (!GUNDAM::_dp_iso::RefineCandidateSet(
                pattern , 
                data_graph, 
                pattern_candidate_set)) {
    return std::make_pair(0, 0);
  }
  return GarSupp<match_semantics>(
                 gar, data_graph, pattern_candidate_set, 
                 x_supp_limit, 
                xy_supp_limit,
                   time_limit,
                   time_limit_per_supp,
                pivoted_vertex_set );
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
inline std::pair<typename DataGraphType::VertexCounterType,
                 typename DataGraphType::VertexCounterType> 
  IncrementalGarSupp(
     GraphAssociationRule<GraphPatternType,
                             DataGraphType>& gar,
                             DataGraphType&  data_graph,
    std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> &delta_target_graph,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
        double     time_limit = -1.0,
        double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using CandidateSet = std::map<GraphPatternVertexHandle, 
                       std::vector<DataGraphVertexHandle>>;
  CandidateSet candidate_set;
  if (!GUNDAM::_dp_iso_using_match::InitCandidateSet<match_semantics>(
           gar.pattern(), 
              data_graph, 
           candidate_set)) {
    return std::pair(0, 0);
  }

  if (!GUNDAM::_dp_iso_using_match::RefineCandidateSet(
           gar.pattern(), 
              data_graph, 
           candidate_set)) {
    return std::pair(0, 0);
  }

  std::pair<typename DataGraphType::VertexCounterType,
            typename DataGraphType::VertexCounterType> gar_match_result;

  std::sort(delta_target_graph.begin(), 
            delta_target_graph.end());

  std::vector<GraphPatternVertexHandle> has_delta_target_graph_pattern_vertex;

  for (auto &[query_handle, target_list] : candidate_set) {
    bool found_new_vertex = false;
    for (auto &target_handle : target_list) {
      if (std::binary_search(delta_target_graph.begin(),
                             delta_target_graph. end (),
                                   target_handle)) {
        found_new_vertex = true;
        break;
      }
    }
    if (!found_new_vertex) {
      continue;
    }
    has_delta_target_graph_pattern_vertex.emplace_back(query_handle);
  }
  std::sort(has_delta_target_graph_pattern_vertex.begin(),
            has_delta_target_graph_pattern_vertex.end());
  int total_mask = (1 << (has_delta_target_graph_pattern_vertex.size()));
  for (int mask = 1; mask < total_mask; mask++) {
    std::vector<GraphPatternVertexHandle> this_mask_vertex;
    for (int bit_pos = 0;
             bit_pos < has_delta_target_graph_pattern_vertex.size(); 
             bit_pos++) {
      if (mask & (1 << bit_pos)) {
        this_mask_vertex.emplace_back(
            has_delta_target_graph_pattern_vertex[bit_pos]);
      }
    }
    CandidateSet copy_candidate_set{candidate_set};
    for (auto &[query_handle, target_list] : copy_candidate_set) {
      if (std::binary_search(this_mask_vertex.begin(), 
                             this_mask_vertex.end(),
                             query_handle)) {
        std::vector<DataGraphVertexHandle> this_vertex_target_list;
        for (auto &target_handle : target_list) {
          if (std::binary_search(delta_target_graph.begin(),
                                 delta_target_graph.end(), 
                                       target_handle)) {
            this_vertex_target_list.emplace_back(target_handle);
          }
        }
        std::swap(this_vertex_target_list, 
                              target_list);
        continue;
      }
      std::vector<DataGraphVertexHandle> this_vertex_target_list;
      for (auto &target_handle : target_list) {
        if (!std::binary_search(delta_target_graph.begin(),
                                delta_target_graph.end(), 
                                      target_handle)) {
          this_vertex_target_list.emplace_back(target_handle);
        }
      }
      std::swap(this_vertex_target_list, 
                            target_list);
    }

    auto partial_result = GarSupp<match_semantics>(
                                  gar, data_graph, copy_candidate_set, 
                                  satisfy_x_supp_callback,
                                 satisfy_xy_supp_callback,
                                    time_limit,
                                    time_limit_per_supp,
                                 pivoted_vertex_set );

    gar_match_result.first  += partial_result.first;
    gar_match_result.second += partial_result.second;
  }
  return gar_match_result;
}

template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
inline std::pair<typename DataGraphType::VertexCounterType,
                 typename DataGraphType::VertexCounterType> 
  IncrementalGarSupp(
     GraphAssociationRule<GraphPatternType,
                             DataGraphType>& gar,
                             DataGraphType&  data_graph,
    std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type> &delta_target_graph,
        int64_t  x_supp_limit = -1,
        int64_t xy_supp_limit = -1,
        double     time_limit = -1.0,
        double     time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using MatchType = GUNDAM::Match<GraphPatternType,
                                     DataGraphType>;

  assert(x_supp_limit == -1 
     || xy_supp_limit == -1 
     || (x_supp_limit >= xy_supp_limit));

  int64_t x_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)>
    satisfy_x_supp_callback 
        = [&x_supp_counter,
           &x_supp_limit](const MatchType& match) {
    if (x_supp_limit == -1){
      // x_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(x_supp_counter < x_supp_limit);
    x_supp_counter++;
    if (x_supp_counter == x_supp_limit){
      // has reached x_supp_limit
      // stop matching
      return false;
    }
    return true;
  };

  int64_t xy_supp_counter = 0;

  std::function<bool(const GUNDAM::Match<GraphPatternType, DataGraphType>&)> 
    satisfy_xy_supp_callback 
        = [&xy_supp_counter,
           &xy_supp_limit](const MatchType& match) {
    if (xy_supp_limit == -1) {
      // xy_supp_limit is not specified 
      // continue matching
      return true;
    }
    assert(xy_supp_counter < xy_supp_limit);
    xy_supp_counter++;
    if (xy_supp_counter == xy_supp_limit) {
      // has reached xy_supp_limit
      // stop matching
      return false;
    }
    return true;
  };
  
  return IncrementalGarSupp(gar, data_graph, delta_target_graph,
                            satisfy_x_supp_callback,
                            satisfy_xy_supp_callback,
                            time_limit,
                            time_limit_per_supp,
                            pivoted_vertex_set);
}


/*  ################################################  *
 *  ##                                            ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##          pivoted_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##      pivoted_rhs_pattern_match_callback    ##  *
 *  ##      // currently this level is omitted    ##  *
 *  ##      // since it must be same with         ##  *
 *  ##      // pivoted_pattern now                ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##common_pivoted_literal_pattern_match_callback#  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##  pivoted_literal_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##              gar_pattern_match_callback    ##  *
 *  ##                                            ##  *
 *  ################################################  */

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<std::vector<typename    DataGraphType::VertexCounterType>,
          std::vector<typename    DataGraphType::VertexCounterType>>
  GarSuppForSet(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& common_gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::vector<gar::LiteralInfo<GraphPatternType,
                         DataGraphType>>& new_literals_to_check,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                                   const bool store_match,
     std::vector<std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>>>& match_vec_for_all,
     std::function<bool(const unsigned, const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const unsigned, const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_xy_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  all_satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> all_satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using VertexIDType = typename GUNDAM::VertexID<GraphPatternType>::type;

  assert(time_limit == -1.0
      || time_limit > 0);

  assert(time_limit_per_supp == -1.0
      || time_limit_per_supp > 0);

  assert(pivoted_vertex_set.empty() 
      || pivoted_vertex_set.size() > gar_pattern_to_data_graph_partial_match.size());

  if (// time_limit is set but time_limit_per_supp is not set
      (time_limit_per_supp == -1.0 && time_limit != -1.0)
      // both time_limit and time_limit_per_supp are set
      // but there are time_limit_per_supp > time_limit
   || (time_limit_per_supp != -1.0 && time_limit != -1.0
    && time_limit_per_supp > time_limit)) {
    time_limit_per_supp = time_limit;
  }
  match_vec_for_all.resize(new_literals_to_check.size());
  auto& common_gar_pattern = common_gar.pattern();

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;


  /* ############################################################## *
   * ##   obtain common pattern contained in all literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
        common_pivoted_literal_gar = _gar_supp::LiteralGar(common_gar, pivoted_vertex_set);
  auto& common_pivoted_literal_pattern = common_pivoted_literal_gar.pattern();

  /* ############################################################## *
   * ##   obtain patterns contained in all literals and pivotes   ## *
   * ############################################################## */
  std::vector<GraphAssociationRule<GraphPatternType, 
                          DataGraphType>>
        gars_vec, pivoted_literal_gars_vec;
  // std::vector<GraphAssociationRule<GraphPatternType, 
  //                         DataGraphType>>
  //       pivoted_literal_gars_to_check_vec;  
  for (unsigned idx = 0; idx < new_literals_to_check.size(); idx++) {
    GraphAssociationRule<GraphPatternType, 
                        DataGraphType>
                    gar_to_check = common_gar;
    gar_to_check.AddX(new_literals_to_check[idx]);
    gars_vec.emplace_back(gar_to_check);
    GraphAssociationRule<GraphPatternType, 
                        DataGraphType>
          pivoted_literal_gar =
                  _gar_supp::LiteralGar(gar_to_check, pivoted_vertex_set);
    pivoted_literal_gars_vec.emplace_back(pivoted_literal_gar);
  }

  /* ############################################################## *
   * ##   obtain pattern contained in rhs literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
         pivoted_rhs_gar = _gar_supp::RhsGar(common_gar, pivoted_vertex_set);
  auto&  pivoted_rhs_pattern = pivoted_rhs_gar.pattern();
  //assert(pivoted_literal_pattern.CountVertex()
  //        >= pivoted_rhs_pattern.CountVertex());
  //assert(GUNDAM::SubGraphOf(pivoted_rhs_pattern,
  //                      pivoted_literal_pattern));

  /* ############################################# *
   * ##   obtain pattern contained in pivotes   ## *
   * ############################################# */
  GraphPatternType pivoted_pattern = pivoted_vertex_set.empty()?
                                     pivoted_rhs_pattern
                                   : GUNDAM::PreserveVertexSet(common_gar.pattern(),
                                     pivoted_vertex_set);
  assert(pivoted_rhs_pattern.CountVertex() 
          >= pivoted_pattern.CountVertex());
  assert(pivoted_vertex_set.empty()
      || pivoted_pattern.CountVertex() 
      == pivoted_vertex_set.size());
  assert(GUNDAM::SubGraphOf(pivoted_pattern,
                        pivoted_rhs_pattern));
  #ifndef NDEBUG
  for (auto map_it = gar_pattern_to_data_graph_partial_match.MapBegin();
           !map_it.IsDone();
            map_it++) {
    assert(pivoted_pattern.FindVertex(map_it->src_handle()->id()));
  }
  #endif // NDEBUG

  /* ########################################## *
   * ##  parameters to be used in callbacks  ## *
   * ########################################## */
  const bool pivoted_pattern_same_as_pivoted_rhs_pattern 
               = GUNDAM::SamePattern(pivoted_pattern,
                                     pivoted_rhs_pattern);
  if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
    std::cout << "pivoted_pattern should be same with pivoted_rhs_pattern" << std::endl;
    return std::pair(std::vector<typename    DataGraphType::VertexCounterType>(),
          std::vector<typename    DataGraphType::VertexCounterType>());
  }

  const bool pivoted_rhs_gar_same_as_common_pivoted_literal_gar
               = SameGar(pivoted_rhs_gar, common_pivoted_literal_gar);
  
  std::vector<unsigned> pivoted_literal_pattern_same_as_common_pivoted_literal_pattern;
  for (unsigned pivoted_literal_gar_idx = 0;
                pivoted_literal_gar_idx < pivoted_literal_gars_vec.size();
                pivoted_literal_gar_idx++) {
    if (GUNDAM::SamePattern(common_pivoted_literal_pattern,
                            pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern())) {
      pivoted_literal_pattern_same_as_common_pivoted_literal_pattern.push_back(
                    pivoted_literal_gar_idx);
    }
  }

  std::vector<bool> pivoted_literal_pattern_same_as_gar_pattern(
          new_literals_to_check.size(), false);
  for (unsigned idx = 0;
          idx < pivoted_literal_pattern_same_as_gar_pattern.size();
          idx++) {
    pivoted_literal_pattern_same_as_gar_pattern[idx]
            = GUNDAM::SamePattern(pivoted_literal_gars_vec[idx].pattern(),
                                  gars_vec[idx].pattern());
  }

  std::set<VertexIDType> common_vertex_id_set;
  for (auto vertex_it = common_pivoted_literal_pattern.VertexBegin();
          !vertex_it.IsDone();
          vertex_it++) {
    common_vertex_id_set.emplace(vertex_it->id());
  }

  std::map<VertexIDType, std::set<VertexIDType>> delta_vid_to_gar_idx;

  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    auto &pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
    for (auto vertex_it = pivoted_literal_pattern.VertexBegin();
              !vertex_it.IsDone();
              vertex_it++) {
      if (common_vertex_id_set.find(vertex_it->id()) == common_vertex_id_set.end()) {
        delta_vid_to_gar_idx[vertex_it->id()].insert(idx);
      }
    }
  }


  std::vector<GraphAssociationRule<GraphPatternType, 
                        DataGraphType>> delta_pivoted_literal_gars_vec;
  std::vector<std::pair<VertexIDType, unsigned>> delta_vid_to_size;
  std::vector<VertexIDType> delta_vid_vec;
  std::vector<bool> visited_gar_idx(pivoted_literal_gars_vec.size(), false);
  std::vector<std::set<unsigned>>
        delta_pivoted_literal_gar_to_pivoted_literal_gar_vec;
  std::vector<std::set<unsigned>>
        delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec;
  std::map<unsigned, unsigned> gar_to_corresponding_delta_pivoted_literal_gar;


  for (auto [vid, idx_set] : delta_vid_to_gar_idx) {
    delta_vid_to_size.emplace_back(vid, idx_set.size());
  }

  // bool delta_pivoted_literal_gar_empty = false;
  // if (delta_vid_to_size.empty()) {
  //   delta_pivoted_literal_gar_empty = true;
  // }

  while (!delta_vid_to_size.empty()) {
    std::sort(delta_vid_to_size.begin(), delta_vid_to_size.end(),
        [](const std::pair<VertexIDType, unsigned> &p1,
           const std::pair<VertexIDType, unsigned> &p2) {
             return p1.second > p2.second; });
    auto front_vid = delta_vid_to_size.front().first;
    //delta_vid_vec.emplace_back(front_id);

    //create delta_pivoted_literal_gar
    auto delta_vid_handle = common_gar_pattern.FindVertex(front_vid);
    if (!delta_vid_handle) {
      std::cout << "should have the handle" << std::endl;
    }
    auto delta_pivoted_vertex_set = pivoted_vertex_set;
    delta_pivoted_vertex_set.emplace_back(delta_vid_handle);
    GraphAssociationRule<GraphPatternType, 
                        DataGraphType>
          delta_pivoted_literal_gar =
                  _gar_supp::LiteralGar(common_gar, delta_pivoted_vertex_set);
    delta_pivoted_literal_gars_vec.emplace_back(delta_pivoted_literal_gar);
    auto &delta_pivoted_literal_pattern = delta_pivoted_literal_gar.pattern();

    std::set<unsigned> delta_pivoted_literal_gar_to_pivoted_literal_gar;
    std::set<unsigned> delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern;
    for (auto idx : delta_vid_to_gar_idx[front_vid]) {
      auto &pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
      if (visited_gar_idx[idx]) {
        continue;
//        std::cout << "a pivoted_literal_gar should not be visited twice" << std::endl;
//        return std::pair(std::vector<uint64_t>(), std::vector<uint64_t>());
      }
      visited_gar_idx[idx] = true;

      // update parental relationship of delta_pivoted_literal gar to pivoted_literal_gar
      if (GUNDAM::SamePattern(pivoted_literal_pattern,
                              delta_pivoted_literal_pattern)) {
        delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern.insert(idx);
      } else {
        delta_pivoted_literal_gar_to_pivoted_literal_gar.insert(idx);
      }

      gar_to_corresponding_delta_pivoted_literal_gar[idx]
              = delta_pivoted_literal_gar_to_pivoted_literal_gar_vec.size();

      // updated the sorted vector delta_vid_to_size
      for (auto vertex_it = pivoted_literal_pattern.VertexBegin();
                !vertex_it.IsDone();
                vertex_it++) {
        if (common_vertex_id_set.find(vertex_it->id()) == common_vertex_id_set.end()) {
          bool found_vid = false;
          for (auto vec_it = delta_vid_to_size.begin();
                    vec_it <  delta_vid_to_size.end();) {
            if (vec_it->first == vertex_it->id()) {
              found_vid = true;
              unsigned current_size = vec_it->second;
              if (current_size == 0) {
                std::cout << "current size should not be 0" << std::endl;
              }
              current_size--;
              if (current_size == 0) {
                vec_it = delta_vid_to_size.erase(vec_it);
                continue;
              }
              *vec_it = std::make_pair(vertex_it->id(), current_size);
              vec_it++;
            } else {
              vec_it++;
            }
            //std::cout << "here 2" << std::endl;
          }
          if (!found_vid) {
            std::cout << "should have found a vid in delta_vid_to_size" << std::endl;
            return std::pair(std::vector<uint64_t>(), std::vector<uint64_t>());
          }
        }
      }
    //std::cout << "here1" << std::endl;
    }
    delta_pivoted_literal_gar_to_pivoted_literal_gar_vec.
          emplace_back(delta_pivoted_literal_gar_to_pivoted_literal_gar);
    delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec.
          emplace_back(delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern);
  }

  unsigned to_check = 0;
  for (unsigned i = 0; i < delta_pivoted_literal_gar_to_pivoted_literal_gar_vec.size(); i++) {
    to_check += delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[i].size();
    to_check += delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[i].size();
  }
  to_check += pivoted_literal_pattern_same_as_common_pivoted_literal_pattern.size();

  if (to_check != new_literals_to_check.size()) {
    std::cout << "size of to check is not same with new literals to check" << std::endl;
  }
    //std::cout << "finish delta tree generation" << std::endl;
  

  /* ############################ *
   * ##  from pivoted_pattern  ## *
   * ############################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match(
    pivoted_rhs_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    common_pivoted_literal_pattern_to_pivoted_pattern_partial_match(
    common_pivoted_literal_pattern,   pivoted_pattern, "same_id_map");

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    delta_pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    delta_pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec.emplace_back(
      delta_pivoted_literal_gars_vec[idx].pattern(), pivoted_pattern, "same_id_map");
  }


  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec.emplace_back(
      pivoted_literal_gars_vec[idx].pattern(), pivoted_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_pivoted_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_pivoted_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, pivoted_pattern, "same_id_map");
  }
  


  /* ################################ *
   * ##  from pivoted_rhs_pattern  ## *
   * ################################ */

  GUNDAM::Match<GraphPatternType, GraphPatternType>
    common_pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec(
      common_pivoted_literal_pattern,   pivoted_rhs_pattern, "same_id_map");

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    delta_pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    delta_pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec.emplace_back(
      delta_pivoted_literal_gars_vec[idx].pattern(), pivoted_rhs_pattern, "same_id_map");
  }  

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec.emplace_back(
      pivoted_literal_gars_vec[idx].pattern(), pivoted_rhs_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_pivoted_rhs_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_pivoted_rhs_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, pivoted_rhs_pattern, "same_id_map");
  }

  /* ########################################### *
   * ##  from common_pivoted_literal_pattern  ## *
   * ########################################### */
  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    delta_pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    delta_pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec.emplace_back(
      delta_pivoted_literal_gars_vec[idx].pattern(), common_pivoted_literal_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec.emplace_back(
      pivoted_literal_gars_vec[idx].pattern(), common_pivoted_literal_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, common_pivoted_literal_pattern, "same_id_map");
  }

  /* ########################################## *
   * ##  from delta_pivoted_literal_pattern  ## *
   * ########################################## */  
  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec(
            pivoted_literal_gars_vec.size());
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto pivoted_literal_gar_idx :
              delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[idx]) {
      pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                pivoted_literal_gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
                  pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern(),
                  delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }

  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto pivoted_literal_gar_idx :
              delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[idx]) {
      pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                pivoted_literal_gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
                  pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern(),
                  delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec(gars_vec.size());
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto gar_idx : // gar_pattern_idx is the same with pivoted_literal_pattern_idx
              delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[idx]) {
      gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
            common_gar_pattern,
            delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }

  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto gar_idx :
              delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[idx]) {
      gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
            common_gar_pattern,
            //gars_vec[gar_idx].pattern(),
            delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }


  /* #################################### *
   * ##  from pivoted_literal_pattern  ## *
   * #################################### */
  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_pivoted_literal_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
      //gars_vec[idx].pattern(), pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
  }


  // GUNDAM::Match<GraphPatternType, GraphPatternType>
  //   gar_pattern_to_pivoted_literal_pattern_partial_match(
  //   gar_pattern,   pivoted_literal_pattern, "same_id_map");


  //####################################
  //* prepare candiate set for new gars*
  //####################################
  CandidateSetContainer pivoted_pattern_candidate_set,
                    pivoted_rhs_pattern_candidate_set,
         common_pivoted_literal_pattern_candidate_set;

  std::vector<CandidateSetContainer> 
          delta_pivoted_literal_pattern_candidate_set(delta_pivoted_literal_gars_vec.size()),
          pivoted_literal_pattern_candidate_set(pivoted_literal_gars_vec.size());

  for (const auto& [gar_pattern_candidate_ptr,
                    gar_pattern_candidate]
                  : gar_pattern_candidate_set_for_vertexes_contained_in_literals) {
    assert(gar_pattern_candidate_ptr);
    if (!gar_pattern_candidate_ptr) {
        std::cout << "error candidate" << std::endl;
    return std::pair(std::vector<typename    DataGraphType::VertexCounterType>(),
          std::vector<typename    DataGraphType::VertexCounterType>());

    }
     auto kCommonPivotedLiteralPatternPtr 
              = common_pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
    if (kCommonPivotedLiteralPatternPtr) {
      // this vertex is contained in common_pivoted_literal_pattern
      auto [ common_pivoted_literal_pattern_candidate_set_it,
             common_pivoted_literal_pattern_candidate_set_ret ]
           = common_pivoted_literal_pattern_candidate_set.emplace(
            kCommonPivotedLiteralPatternPtr,
            gar_pattern_candidate);
      // should added successfully
      assert(common_pivoted_literal_pattern_candidate_set_ret);

      for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
        auto& delta_pivoted_literal_pattern = delta_pivoted_literal_gars_vec[idx].pattern();
        auto kDeltaPivotedLiteralPatternPtr
              = delta_pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        auto [delta_pivoted_literal_pattern_candidate_set_it,
              delta_pivoted_literal_pattern_candidate_set_ret ]
            = delta_pivoted_literal_pattern_candidate_set[idx].emplace(
              kDeltaPivotedLiteralPatternPtr,
              gar_pattern_candidate);
        assert(delta_pivoted_literal_pattern_candidate_set_ret);
        // should added successfully
      }

      for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
        auto& pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
        auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        auto [ pivoted_literal_pattern_candidate_set_it,
               pivoted_literal_pattern_candidate_set_ret ]
            = pivoted_literal_pattern_candidate_set[idx].emplace(
              kPivotedLiteralPatternPtr,
              gar_pattern_candidate);
        assert(pivoted_literal_pattern_candidate_set_ret);
        // should added successfully
      }

      const auto kPivotedRhsPatternPtr
                = pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id());
      if (kPivotedRhsPatternPtr) {
        // this vertex is contained in pivoted_rhs_pattern
        auto [ pivoted_rhs_pattern_candidate_set_it,
               pivoted_rhs_pattern_candidate_set_ret ]
             = pivoted_rhs_pattern_candidate_set.emplace(
              kPivotedRhsPatternPtr,
                gar_pattern_candidate);
        // should added successfully
        assert(pivoted_rhs_pattern_candidate_set_ret);
        const auto kPivotedPatternPtr
                  = pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedPatternPtr){
          // this vertex is contained in pivoted_rhs_pattern
          auto [ pivoted_pattern_candidate_set_it,
                 pivoted_pattern_candidate_set_ret ]
               = pivoted_pattern_candidate_set.emplace(
                 kPivotedPatternPtr,
                  gar_pattern_candidate);
          // should added successfully
          assert(pivoted_pattern_candidate_set_ret);
        }
      }
    } else {
      for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
        auto& delta_pivoted_literal_pattern = delta_pivoted_literal_gars_vec[idx].pattern();
        auto kDeltaPivotedLiteralPatternPtr
              = delta_pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kDeltaPivotedLiteralPatternPtr) {
          auto [ delta_pivoted_literal_pattern_candidate_set_it,
                delta_pivoted_literal_pattern_candidate_set_ret ]
              = delta_pivoted_literal_pattern_candidate_set[idx].emplace(
                kDeltaPivotedLiteralPatternPtr,
                gar_pattern_candidate);
          // should added successfully
          assert(delta_pivoted_literal_pattern_candidate_set_ret);
        }
      }

      for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
        auto& pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
        auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedLiteralPatternPtr) {
          auto [ pivoted_literal_pattern_candidate_set_it,
                pivoted_literal_pattern_candidate_set_ret ]
              = pivoted_literal_pattern_candidate_set[idx].emplace(
                kPivotedLiteralPatternPtr,
                gar_pattern_candidate);
          // should added successfully
          assert(pivoted_literal_pattern_candidate_set_ret);
        }
      }
    }
  }

  if (new_literals_to_check.size() == 0) {
    std::cout << "literal size is 0" << std::endl;
  }

  std::vector<bool> satisfy_x_callback_gar(gars_vec.size(), false),
                    satisfy_xy_callback_gar(gars_vec.size(), false),
                    satisfy_x_callback_delta_gar(
                        delta_pivoted_literal_gars_vec.size(), false);
  
  std::vector<typename DataGraphType::VertexCounterType>
                        x_supp_counter_vec(gars_vec.size(), 0),
                       xy_supp_counter_vec(gars_vec.size(), 0);

  std::vector<unsigned> satisfy_x_callback_delta_gar_num(
                              delta_pivoted_literal_gars_vec.size(), 0);

  std::vector<bool> has_x_match_gar(gars_vec.size(), false),
                    all_has_x_match_delta_pivoted_literal_gar(
                          delta_pivoted_literal_gars_vec.size(), false);
  
  // std::vector<bool> gar_satisfy_x_supp(gars_vec.size(), false),
  //                   gar_has_match_vec(gars_vec.size(), false);
  //                   delta_pivoted_literal_gar_has_match(
  //                         delta_pivoted_literal_gars_vec.size(), false),
                    
  unsigned gar_idx = 0,
           pivoted_literal_gar_idx = 0,
           delta_pivoted_literal_gar_idx = 0;
  
  bool all_has_satisfy_x = false;

  // bool has_satisfy_x       = false,
  //      has_satisfy_x_not_y = false,
  //      has_satisfy_x_and_y = false;

  // typename DataGraphType::VertexCounterType x_supp = 0,
  //                                          xy_supp = 0;

  /* ################################################ *
   * ##  callbacks from gar_pattern to data graph  ## *
   * ################################################ */
  std::function<bool(const MatchType&)> 
    gar_pattern_prune_callback
    = [](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // prune nothing, continue the matching
    return false;
  };

  // based on the partial match from literal_kPattern to data graph
  // complete the match of the entire gar_pattern 
  std::function<bool(const MatchType&)> 
    gar_pattern_match_callback 
    = [&all_has_satisfy_x,
       &has_x_match_gar,
       &gar_idx,
       &x_supp_counter_vec,
       &store_match,
       &match_vec_for_all](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // only need to complete the match, find one match that 
    // satisfy all literal is enough
    all_has_satisfy_x = true;
    x_supp_counter_vec[gar_idx]++;
    has_x_match_gar[gar_idx] = true;

    if (store_match) {
      match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_match);
    }

    // terminate matching
    return false;
  };

  /* ############################################################ *
   * ##  callbacks from pivoted_literal_pattern to data graph  ## *
   * ############################################################ */
  // if the current pivoted_literal_pattern_to_data_graph_match
  // has already violate a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr :
            pivoted_literal_gars_vec[pivoted_literal_gar_idx].x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_match)){
        // this literal is not considered in the current
        // pivoted_literal_pattern_to_data_graph_match,
        // this match cannot be pruned
        //
        // move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        //std::cout << "exist x not match pivoted" << std::endl;
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };


  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    all_has_satisfy_x = false;
    gar_idx = pivoted_literal_gar_idx;
    assert(!has_x_match_gar[gar_idx]);

    if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
      // pivoted_rhs_gar is not the same as pivoted_literal_gar
      // first match the pivoted_rhs_literal, then match
      // the pivoted_literal_gar

      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
        = pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_pivoted_literal_pattern_partial_match_vec[
                            gar_idx]);
      assert(!all_has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              common_gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
      if (all_has_satisfy_x) {
        assert(has_x_match_gar[gar_idx]);
        //has_x_match_gar[gar_idx] = true;
        // has found satisfy X, terminate matching
        return false;
      }
      // continue matching
      return true;
    }
    assert(!all_has_satisfy_x);
    // gar_pattern of literal_gar is the same as the
    // entire gar_pattern, does not need further match
    all_has_satisfy_x = true;

    if (store_match) {
      GUNDAM::Match<GraphPatternType, DataGraphType>
                      gar_pattern_to_data_graph_entire_match 
        = pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_pivoted_literal_pattern_partial_match_vec[
                            gar_idx]);
      match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
    }

    x_supp_counter_vec[gar_idx]++;
    has_x_match_gar[gar_idx] = true;
    // has found satisfy X, terminate matching
    return false;
  };

  std::function<bool(const MatchType&)>
    delta_pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    delta_pivoted_literal_pattern_to_data_graph_match) -> bool {
    //prune nothing
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    delta_pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    delta_pivoted_literal_pattern_to_data_graph_match) -> bool {
    unsigned gar_num_no_need_check = 0;

    if (all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]) {
      std::cout << "should not be here delta" << std::endl;
    }

    assert(!all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]);

    for (auto current_pivoted_literal_gar_idx :
            delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[
                  delta_pivoted_literal_gar_idx]) {
      pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
      gar_idx = pivoted_literal_gar_idx;

      if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
        gar_num_no_need_check++;
        continue;
      }
      all_has_satisfy_x = false;
      
      if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
                  delta_pivoted_literal_pattern_to_data_graph_match)) {
        continue;
      }

      if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
        GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
          = delta_pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                          gar_idx]);

        GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              common_gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
        if (all_has_satisfy_x) {
          gar_num_no_need_check++;
        }
      } else {
        // continue matching
        if (store_match) {
          GUNDAM::Match<GraphPatternType, DataGraphType>
                          gar_pattern_to_data_graph_entire_match 
            = delta_pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                          gar_idx]);
          match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
        }

        x_supp_counter_vec[gar_idx]++;
        has_x_match_gar[gar_idx] = true;
        gar_num_no_need_check++;
      }
    }

    for (auto current_pivoted_literal_gar_idx :
              delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[
                  delta_pivoted_literal_gar_idx]) {
      pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
      gar_idx = pivoted_literal_gar_idx;
      if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
        gar_num_no_need_check++;
        continue;
      }

      all_has_satisfy_x = false;
      auto& pivoted_literal_pattern
              = pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern();
      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      pivoted_literal_pattern_to_data_graph_partial_match
        = delta_pivoted_literal_pattern_to_data_graph_match(
                pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                    pivoted_literal_gar_idx]);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph,
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set[pivoted_literal_gar_idx],
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
      if (all_has_satisfy_x) {
        gar_num_no_need_check++;
      }
    }

    all_has_satisfy_x = false;
    if (gar_num_no_need_check
                  == (delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[
                      delta_pivoted_literal_gar_idx].size()
                  + delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[
                      delta_pivoted_literal_gar_idx].size())) {
      all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx] = true;
      all_has_satisfy_x = true;
      return false;
    }
    return true;
  };


  /* ############################################################ *
   * ##  callbacks from common_pivoted_literal_pattern to data graph  ## *
   * ############################################################ */
  // if the current pivoted_literal_pattern_to_data_graph_match
  // has already violate a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    common_pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    common_pivoted_literal_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : common_pivoted_literal_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(common_pivoted_literal_pattern_to_data_graph_match)){
        // this literal is not considered in the current
        // pivoted_literal_pattern_to_data_graph_match,
        // this match cannot be pruned
        //
        // move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(common_pivoted_literal_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        //std::cout << "exist x not match common" << std::endl;
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };


  std::function<bool(const MatchType&)>
    common_pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    common_pivoted_literal_pattern_to_data_graph_match) -> bool {
    unsigned gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback = 0;
    for (auto current_pivoted_literal_gar_idx :
            pivoted_literal_pattern_same_as_common_pivoted_literal_pattern) {
      pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
      gar_idx = pivoted_literal_gar_idx;

      if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
        continue;
      }

      all_has_satisfy_x = false;

      if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
            common_pivoted_literal_pattern_to_data_graph_match)) {
        continue;
      }

      if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
        GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
          = common_pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec[
                          gar_idx]);

        GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              common_gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
        if (all_has_satisfy_x) {
          gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
        }
      } else {
        if (store_match) {
          GUNDAM::Match<GraphPatternType, DataGraphType>
                          gar_pattern_to_data_graph_entire_match 
            = common_pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec[
                          gar_idx]);
          match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
        }

        has_x_match_gar[gar_idx] = true;
        x_supp_counter_vec[gar_idx]++;
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
      }
    }

    for (delta_pivoted_literal_gar_idx = 0;
         delta_pivoted_literal_gar_idx < delta_pivoted_literal_gars_vec.size();
         delta_pivoted_literal_gar_idx++) {
    // delta_pivoted_literal_gar can not be the same with common_pivoted_literal_gar

      if (satisfy_x_callback_delta_gar[delta_pivoted_literal_gar_idx]
              || all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]) {
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
        continue;
      }
      auto& delta_pivoted_literal_pattern
              = delta_pivoted_literal_gars_vec[delta_pivoted_literal_gar_idx].pattern();
      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      delta_pivoted_literal_pattern_to_data_graph_partial_match
        = common_pivoted_literal_pattern_to_data_graph_match(
                delta_pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec[
                    delta_pivoted_literal_gar_idx]);
      all_has_satisfy_x = false;
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              delta_pivoted_literal_pattern,   data_graph,
                              delta_pivoted_literal_pattern_to_data_graph_partial_match,
                              delta_pivoted_literal_pattern_candidate_set[delta_pivoted_literal_gar_idx],
                              delta_pivoted_literal_pattern_prune_callback,
                              delta_pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
      if (all_has_satisfy_x) {
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
      }
    }
    
    all_has_satisfy_x = false;
    if (gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback
            == delta_pivoted_literal_gars_vec.size()
              + pivoted_literal_pattern_same_as_common_pivoted_literal_pattern.size()) {
      all_has_satisfy_x = true;
      return false;
    }
    return true;
  };



  /* ######################################################## *
   * ##  callbacks from pivoted_rhs_pattern to data graph  ## *
   * ######################################################## */
  // if the current match_state has already violate
  // a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_prune_callback
    = [&](const MatchType&
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    std::cout << "should not use this prune callack" << std::endl;
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_match_callback
    = [&](const MatchType& 
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    std::cout << "should not use this match callback" << std::endl;
    return true;
  };

  /* #################################################### *
   * ##  callbacks from pivoted_pattern to data graph  ## *
   * #################################################### */
  std::function<bool(const MatchType&)>
    pivoted_pattern_prune_callback
    = [](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    // prune nothing
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_pattern_match_callback
    = [&](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    /* ######################################## *
     * ##  get a match for pivot vertex set  ## *
     * ######################################## */
    all_has_satisfy_x       = false;
    has_x_match_gar.clear();
    has_x_match_gar.resize(gars_vec.size(), false);
    all_has_x_match_delta_pivoted_literal_gar.clear();
    all_has_x_match_delta_pivoted_literal_gar.resize(
              delta_pivoted_literal_gars_vec.size(), false);
    
    /* ####################################################### *
     * ##  get the paritial match from pivoted_rhs_pattern  ## *
     * ##  to data_graph                                    ## *
     * ####################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_rhs_pattern_to_data_graph_partial_match
        = pivoted_pattern_to_data_graph_match(
          pivoted_rhs_pattern_to_pivoted_pattern_partial_match);

    if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar
      std::cout << "it has not been implemented and should not be here" << std::endl;
      return true;
    }
    /* ########################################################### *
     * ##  pivoted_pattern is the same as pivoted_rhs_pattern,  ## *
     * ##  directly whether the match satisfy y literals here   ## *
     * ########################################################### */
    assert(pivoted_rhs_pattern_to_data_graph_partial_match.size()
            == pivoted_pattern_to_data_graph_match.size());

    /* ########################################################### *
     * ##  get the paritial match from pivoted_literal_pattern  ## *
     * ##  to data_graph                                        ## *
     * ########################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    common_pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
            common_pivoted_literal_pattern_to_pivoted_pattern_partial_match);


    if (!pivoted_rhs_gar_same_as_common_pivoted_literal_gar) {
      //std::cout << "rhs not same as common" << std::endl;
      assert( GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_pattern));
      assert(!        SameGar    (pivoted_rhs_gar,     common_pivoted_literal_gar));
      assert(!GUNDAM::SamePattern(pivoted_rhs_pattern, common_pivoted_literal_pattern));

      assert(!all_has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              common_pivoted_literal_pattern,   data_graph, 
                              common_pivoted_literal_pattern_to_data_graph_partial_match,
                              common_pivoted_literal_pattern_candidate_set,
                              common_pivoted_literal_pattern_prune_callback,
                              common_pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    } else {
      //pivoted_rhs_gar is the same as common_pivoted_literal_gar
      assert(        SameGar    (pivoted_rhs_gar, common_pivoted_literal_gar));
      assert(GUNDAM::SamePattern(pivoted_pattern, common_pivoted_literal_pattern));

      // verify x_literals for common_pivoted_literal_gar
      bool satisify_x_literals = true;
      for (const auto& x_literal_ptr : common_pivoted_literal_gar.x_literal_set()) {
        assert(x_literal_ptr->MappedBy(common_pivoted_literal_pattern_to_data_graph_partial_match));
        if (x_literal_ptr->Satisfy(common_pivoted_literal_pattern_to_data_graph_partial_match)) {
          continue;
        }
        satisify_x_literals = false;
        break;
      }
      if (!satisify_x_literals) {
        // continue matching
        //std::cout << "exist x not match" << std::endl;
        assert(!satisify_x_literals);
        return true;
      }


      for (auto current_pivoted_literal_gar_idx :
              pivoted_literal_pattern_same_as_common_pivoted_literal_pattern) {
        pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
        gar_idx = pivoted_literal_gar_idx;

        if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
          continue;
        }

        // GUNDAM::Match<GraphPatternType, DataGraphType>  
        //                 pivoted_literal_pattern_to_data_graph_partial_match
        //   = pivoted_pattern_to_data_graph_match(
        //             pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec[
        //                 pivoted_literal_gar_idx]);

        // if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
        //           pivoted_literal_pattern_to_data_graph_partial_match)) {
        //   continue;
        // }


        if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
                  pivoted_pattern_to_data_graph_match)) {
          continue;
        }

        if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
          GUNDAM::Match<GraphPatternType, DataGraphType>  
                        gar_pattern_to_data_graph_partial_match
            = pivoted_pattern_to_data_graph_match(
                        gar_pattern_to_pivoted_pattern_partial_match_vec[
                            gar_idx]);

          GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kMerge>(
                                common_gar_pattern,   data_graph,
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
        } else {
          if (store_match) {
            GUNDAM::Match<GraphPatternType, DataGraphType>
                            gar_pattern_to_data_graph_entire_match 
              = pivoted_pattern_to_data_graph_match(
                        gar_pattern_to_pivoted_pattern_partial_match_vec[
                            gar_idx]);
            match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
          }

          has_x_match_gar[gar_idx] = true;
          x_supp_counter_vec[gar_idx]++;
        }
      }


      for (delta_pivoted_literal_gar_idx = 0;
          delta_pivoted_literal_gar_idx < delta_pivoted_literal_gars_vec.size();
          delta_pivoted_literal_gar_idx++) {
      // delta_pivoted_literal_gar can not be the same with common_pivoted_literal_gar

        if (satisfy_x_callback_delta_gar[delta_pivoted_literal_gar_idx]
              || all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]) {
          continue;
        }
        auto& delta_pivoted_literal_pattern
                = delta_pivoted_literal_gars_vec[delta_pivoted_literal_gar_idx].pattern();
        GUNDAM::Match<GraphPatternType, DataGraphType>  
                        delta_pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
                  delta_pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec[
                      delta_pivoted_literal_gar_idx]);
        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                delta_pivoted_literal_pattern,   data_graph,
                                delta_pivoted_literal_pattern_to_data_graph_partial_match,
                                delta_pivoted_literal_pattern_candidate_set[delta_pivoted_literal_gar_idx],
                                delta_pivoted_literal_pattern_prune_callback,
                                delta_pivoted_literal_pattern_match_callback,
                                time_limit_per_supp);
      }
    }

     GUNDAM::Match<GraphPatternType, DataGraphType>
           gar_pattern_to_data_graph_partial_match
     = pivoted_pattern_to_data_graph_match(
           gar_pattern_to_pivoted_pattern_partial_match_vec[0]);


    bool exist_x_satisfy = false;
    for (unsigned gars_vec_idx = 0;
                  gars_vec_idx < gars_vec.size();
                  gars_vec_idx++) {
      if (has_x_match_gar[gars_vec_idx]) {
        exist_x_satisfy = true;

        const bool kSatisfyXSuppCallbackRet
              = satisfy_x_supp_callback(gars_vec_idx, gar_pattern_to_data_graph_partial_match);
        if (!kSatisfyXSuppCallbackRet) {
          satisfy_x_callback_gar[gars_vec_idx] = true;
        }
      }
    }

    if (!exist_x_satisfy) {
      return true;
    }

    /* ######################################################## *
    * ##  chech whether the match satisfy rhs literal here  ## *
    * ######################################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_partial_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_partial_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }
    const bool kAllSatisfyXSuppCallbackRet
          = all_satisfy_x_supp_callback(gar_pattern_to_data_graph_partial_match);
    if (!satisify_y_literals) {
      return kAllSatisfyXSuppCallbackRet;
    }
    // satisfy y literals

    for (unsigned gars_vec_idx = 0;
                  gars_vec_idx < gars_vec.size();
                  gars_vec_idx++) {
      if (has_x_match_gar[gars_vec_idx]) {
        xy_supp_counter_vec[gars_vec_idx]++;
        const bool kSatisfyXYSuppCallbackRet
              = satisfy_xy_supp_callback(gars_vec_idx, gar_pattern_to_data_graph_partial_match);
        if (!kSatisfyXYSuppCallbackRet) {
          satisfy_x_callback_gar[gars_vec_idx] = true;
          if (gar_to_corresponding_delta_pivoted_literal_gar.find(gars_vec_idx)
                  == gar_to_corresponding_delta_pivoted_literal_gar.end()) {
            continue;
          }
          unsigned corresponding_delta_pivoted_literal_gar_idx
                = gar_to_corresponding_delta_pivoted_literal_gar[gars_vec_idx];
          satisfy_x_callback_delta_gar_num[corresponding_delta_pivoted_literal_gar_idx]++;
          if (satisfy_x_callback_delta_gar_num[corresponding_delta_pivoted_literal_gar_idx]
                < delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[
                      corresponding_delta_pivoted_literal_gar_idx].size()
                  + delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[
                      corresponding_delta_pivoted_literal_gar_idx].size()) {
            continue;
          }
          satisfy_x_callback_delta_gar[corresponding_delta_pivoted_literal_gar_idx] = true;
        }
      }
    }
    // continue matching
    const bool kAllSatisfyXYSuppCallbackRet
              = all_satisfy_xy_supp_callback(gar_pattern_to_data_graph_partial_match);
    return kAllSatisfyXYSuppCallbackRet
            && kAllSatisfyXSuppCallbackRet;

  };

  MatchType pivoted_pattern_partial_match;
  for (auto pivoted_pattern_vertex_it = pivoted_pattern.VertexBegin();
           !pivoted_pattern_vertex_it.IsDone();
            pivoted_pattern_vertex_it++) {
    auto target_handle = gar_pattern_to_data_graph_partial_match.MapTo(pivoted_pattern_vertex_it->id());
    if (!target_handle) {
      continue;
    }
    pivoted_pattern_partial_match.AddMap(pivoted_pattern_vertex_it,
                                         target_handle);
  }

  // first match pivoted gar_pattern, then match rhs gar_pattern
  GUNDAM::MatchUsingMatch<match_semantics,
                          GUNDAM::MatchAlgorithm::kDagDp,
                          GUNDAM::MergeNecConfig::kNotMerge>(
                          pivoted_pattern, data_graph,
                          pivoted_pattern_partial_match,
                          pivoted_pattern_candidate_set,
                          pivoted_pattern_prune_callback,
                          pivoted_pattern_match_callback,
                          time_limit);
  //assert(x_supp >= xy_supp);
  return std::pair(x_supp_counter_vec, xy_supp_counter_vec);
}


// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  IncGarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
     std::map<typename GUNDAM::VertexHandle<DataGraphType>::type,
              std::set<typename GUNDAM::VertexHandle<DataGraphType>::type>>& support_candidate_set,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  assert(time_limit == -1.0
      || time_limit > 0);

  assert(time_limit_per_supp == -1.0
      || time_limit_per_supp > 0);

  assert(pivoted_vertex_set.empty() 
      || pivoted_vertex_set.size() > gar_pattern_to_data_graph_partial_match.size());

  if (// time_limit is set but time_limit_per_supp is not set
      (time_limit_per_supp == -1.0 && time_limit != -1.0)
      // both time_limit and time_limit_per_supp are set
      // but there are time_limit_per_supp > time_limit
   || (time_limit_per_supp != -1.0 && time_limit != -1.0
    && time_limit_per_supp > time_limit)) {
    time_limit_per_supp = time_limit;
  }

  auto& gar_pattern = gar.pattern();

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;

  /* ############################################################## *
   * ##   obtain pattern contained in all literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
        pivoted_literal_gar = _gar_supp::LiteralGar(gar, pivoted_vertex_set);
  auto& pivoted_literal_pattern = pivoted_literal_gar.pattern();

  /* ############################################################## *
   * ##   obtain pattern contained in rhs literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
         pivoted_rhs_gar = _gar_supp::RhsGar(gar, pivoted_vertex_set);
  auto&  pivoted_rhs_pattern = pivoted_rhs_gar.pattern();
  assert(pivoted_literal_pattern.CountVertex()
          >= pivoted_rhs_pattern.CountVertex());
  assert(GUNDAM::SubGraphOf(pivoted_rhs_pattern,
                        pivoted_literal_pattern));

  /* ############################################# *
   * ##   obtain pattern contained in pivotes   ## *
   * ############################################# */
  GraphPatternType pivoted_pattern = pivoted_vertex_set.empty()?
                                     pivoted_rhs_pattern
                                   : GUNDAM::PreserveVertexSet(gar.pattern(),
                                     pivoted_vertex_set);
  assert(pivoted_rhs_pattern.CountVertex() 
          >= pivoted_pattern.CountVertex());
  assert(pivoted_vertex_set.empty()
      || pivoted_pattern.CountVertex() 
      == pivoted_vertex_set.size());
  assert(GUNDAM::SubGraphOf(pivoted_pattern,
                        pivoted_rhs_pattern));
  #ifndef NDEBUG
  for (auto map_it = gar_pattern_to_data_graph_partial_match.MapBegin();
           !map_it.IsDone();
            map_it++) {
    assert(pivoted_pattern.FindVertex(map_it->src_handle()->id()));
  }
  #endif // NDEBUG

  /* ########################################## *
   * ##  parameters to be used in callbacks  ## *
   * ########################################## */
  const bool pivoted_pattern_same_as_pivoted_rhs_pattern 
               = GUNDAM::SamePattern(pivoted_pattern,
                                     pivoted_rhs_pattern);

  const bool pivoted_rhs_gar_same_as_pivoted_literal_gar 
                           = SameGar(pivoted_rhs_gar,
                                     pivoted_literal_gar);

  const bool pivoted_literal_pattern_same_as_gar_pattern
           = GUNDAM::SamePattern(pivoted_literal_pattern,
                                             gar_pattern);

  /* ############################ *
   * ##  from pivoted_pattern  ## *
   * ############################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match(
    pivoted_rhs_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_pattern_partial_match(
    gar_pattern,   pivoted_pattern, "same_id_map");

  /* ################################ *
   * ##  from pivoted_rhs_pattern  ## *
   * ################################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_rhs_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_rhs_pattern_partial_match(
    gar_pattern,   pivoted_rhs_pattern, "same_id_map");

  /* #################################### *
   * ##  from pivoted_literal_pattern  ## *
   * #################################### */
  GUNDAM::Match<GraphPatternType, GraphPatternType>
    gar_pattern_to_pivoted_literal_pattern_partial_match(
    gar_pattern,   pivoted_literal_pattern, "same_id_map");

  CandidateSetContainer pivoted_pattern_candidate_set,
                    pivoted_rhs_pattern_candidate_set,
                pivoted_literal_pattern_candidate_set;

  // process pivoted_literal_pattern_candidate_set
  //         and pivoted_rhs_pattern_candidate_set
  for (const auto& [gar_pattern_candidate_ptr,
                    gar_pattern_candidate]
                  : gar_pattern_candidate_set_for_vertexes_contained_in_literals) {
    // should not be null
    assert(gar_pattern_candidate_ptr);
    const auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
    if (kPivotedLiteralPatternPtr) {
      // this vertex is contained in pivoted_literal_pattern
      auto [ pivoted_literal_pattern_candidate_set_it,
             pivoted_literal_pattern_candidate_set_ret ]
           = pivoted_literal_pattern_candidate_set.emplace(
            kPivotedLiteralPatternPtr,
            gar_pattern_candidate);
      // should added successfully
      assert(pivoted_literal_pattern_candidate_set_ret);
      const auto kPivotedRhsPatternPtr
                = pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id());
      if (kPivotedRhsPatternPtr){
        // this vertex is contained in pivoted_rhs_pattern
        auto [ pivoted_rhs_pattern_candidate_set_it,
               pivoted_rhs_pattern_candidate_set_ret ]
             = pivoted_rhs_pattern_candidate_set.emplace(
              kPivotedRhsPatternPtr,
                gar_pattern_candidate);
        // should added successfully
        assert(pivoted_rhs_pattern_candidate_set_ret);
        const auto kPivotedPatternPtr
                  = pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedPatternPtr){
          // this vertex is contained in pivoted_rhs_pattern
          auto [ pivoted_pattern_candidate_set_it,
                 pivoted_pattern_candidate_set_ret ]
               = pivoted_pattern_candidate_set.emplace(
                 kPivotedPatternPtr,
                  gar_pattern_candidate);
          // should added successfully
          assert(pivoted_pattern_candidate_set_ret);
        }
      }
      #ifndef NDEBUG
      else {
        // should not be contained in pivoted_rhs_pattern
        assert(!pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      }
      #endif // NDEBUG
    }
    #ifndef NDEBUG
    else {
      // should not be contained in pivoted_rhs_pattern and pivoted_pattern
      assert(!pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      assert(    !pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
    }
    #endif // NDEBUG
  }

  bool has_satisfy_x       = false,
       has_satisfy_x_not_y = false,
       has_satisfy_x_and_y = false;

  typename DataGraphType::VertexCounterType x_supp = 0,
                                           xy_supp = 0;

  const auto kRhsLiteralInfo
         = (*gar.y_literal_set()
                .begin())->info();

  auto rhs_literal_x_id = kRhsLiteralInfo.x_id();
  auto rhs_literal_y_id = kRhsLiteralInfo.y_id();


  /* ################################################ *
   * ##  callbacks from gar_pattern to data graph  ## *
   * ################################################ */
  std::function<bool(const MatchType&)> 
    gar_pattern_prune_callback
    = [](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // prune nothing, continue the matching
    return false;
  };
  // based on the partial match from literal_kPattern to data graph
  // complete the match of the entire gar_pattern 
  std::function<bool(const MatchType&)> 
    gar_pattern_match_callback 
    = [&has_satisfy_x,
       &rhs_literal_x_id,
       &rhs_literal_y_id,
       &support_candidate_set](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // only need to complete the match, find one match that 
    // satisfy all literal is enough
    has_satisfy_x = true;
    auto data_graph_x_handle = gar_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
    auto data_graph_y_handle = gar_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

    support_candidate_set[data_graph_x_handle].insert(data_graph_y_handle);
    // terminate matching
    return false;
  };

  /* ############################################################ *
   * ##  callbacks from pivoted_literal_pattern to data graph  ## *
   * ############################################################ */
  // if the current pivoted_literal_pattern_to_data_graph_match
  // has already violate a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : pivoted_literal_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_match)){
        // this literal is not considered in the current
        // pivoted_literal_pattern_to_data_graph_match,
        // this match cannot be pruned
        //
        // move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    if (!pivoted_literal_pattern_same_as_gar_pattern) {
      // pivoted_rhs_gar is not the same as pivoted_literal_gar
      // first match the pivoted_rhs_literal, then match
      // the pivoted_literal_gar
      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
        = pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_pivoted_literal_pattern_partial_match);
      assert(!has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
      if (has_satisfy_x) {
        // has found satisfy X, terminate matching
        return false;
      }
      // continue matching
      return true;
    }
    assert(!has_satisfy_x);
    // gar_pattern of literal_gar is the same as the
    // entire gar_pattern, does not need further match
    auto data_graph_x_handle = pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
    auto data_graph_y_handle = pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

    support_candidate_set[data_graph_x_handle].insert(data_graph_y_handle);
    has_satisfy_x = true;
    // has found satisfy X, terminate matching
    return false;
  };

  /* ######################################################## *
   * ##  callbacks from pivoted_rhs_pattern to data graph  ## *
   * ######################################################## */
  // if the current match_state has already violate
  // a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_prune_callback
    = [&](const MatchType&
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : pivoted_rhs_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match)) {
        // this literal is not considered in the current
        // pivoted_rhs_pattern_to_data_graph_match, 
        // this match cannot be pruned move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_match_callback
    = [&](const MatchType& 
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    assert(!has_satisfy_x);
    assert(!has_satisfy_x_not_y);
    /* ############################################## *
     * ##  first check whether satisfy y literals  ## *
     * ############################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }

    if (!pivoted_rhs_gar_same_as_pivoted_literal_gar) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar

      /* ####################################################### *
       * ##  get the paritial match from pivoted_rhs_pattern  ## *
       * ##  to data_graph                                    ## *
       * ####################################################### */
      GUNDAM::Match<GraphPatternType, DataGraphType>
      pivoted_literal_pattern_to_data_graph_partial_match
        = pivoted_rhs_pattern_to_data_graph_match(
      pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match);
      
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph, 
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set,
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    }
    else {
      // the pivoted_rhs_gar is the same as pivoted_literal_gar
      // directly match literal_gar to data_graph
      if (!pivoted_literal_pattern_same_as_gar_pattern) {
        /* ####################################################### *
         * ##  get the paritial match from pivoted_rhs_pattern  ## *
         * ##  to data_graph                                    ## *
         * ####################################################### */
        GUNDAM::Match<GraphPatternType, DataGraphType>
                  gar_pattern_to_data_graph_partial_match
        = pivoted_rhs_pattern_to_data_graph_match(
                  gar_pattern_to_pivoted_rhs_pattern_partial_match);

        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                gar_pattern,   data_graph, 
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
      }
      else {
        // pivoted_rhs_gar_pattern is the same as gar_pattern
      auto data_graph_x_handle = pivoted_rhs_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
      auto data_graph_y_handle = pivoted_rhs_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

      support_candidate_set[data_graph_x_handle].insert(data_graph_y_handle);
        has_satisfy_x = true;
      }
    }
    assert(!has_satisfy_x_not_y);
    assert(!has_satisfy_x_and_y);
    if (!has_satisfy_x) {
      // continue matching
      return true;
    }
    // satisfy x
    // reset parameter
    has_satisfy_x = false;
    if (satisify_y_literals) {
      // satisfy x and y
      assert(!has_satisfy_x_not_y);
      has_satisfy_x_and_y = true;
      // continue matching
      return true;
    }
    // satisfy x but not y
    has_satisfy_x_and_y = false;
    has_satisfy_x_not_y =  true;
    // no longer needs matching
    return false;
  };

  /* #################################################### *
   * ##  callbacks from pivoted_pattern to data graph  ## *
   * #################################################### */
  std::function<bool(const MatchType&)>
    pivoted_pattern_prune_callback
    = [&rhs_literal_x_id,
       &rhs_literal_y_id,
       &support_candidate_set](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    auto data_graph_x_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
    
    auto data_graph_y_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);

    if ((!data_graph_y_handle) || (!data_graph_x_handle)) {
      return false;
    }
    if ((support_candidate_set.find(data_graph_x_handle) != support_candidate_set.end())
        && (support_candidate_set[data_graph_x_handle].find(data_graph_y_handle)
              != support_candidate_set[data_graph_x_handle].end())) {
      // the pivot pair has already been verified as legal
      return true;
    }
    // prune nothing
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_pattern_match_callback
    = [&](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    /* ######################################## *
     * ##  get a match for pivot vertex set  ## *
     * ######################################## */
    has_satisfy_x       = false;
    has_satisfy_x_not_y = false;
    has_satisfy_x_and_y = false;
    
    /* ####################################################### *
     * ##  get the paritial match from pivoted_rhs_pattern  ## *
     * ##  to data_graph                                    ## *
     * ####################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_rhs_pattern_to_data_graph_partial_match
      = pivoted_pattern_to_data_graph_match(
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match);

    if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_rhs_pattern,   data_graph, 
                              pivoted_rhs_pattern_to_data_graph_partial_match,
                              pivoted_rhs_pattern_candidate_set,
                              pivoted_rhs_pattern_prune_callback,
                              pivoted_rhs_pattern_match_callback,
                              time_limit_per_supp);
      assert(!has_satisfy_x);
      assert(!has_satisfy_x_and_y 
          || !has_satisfy_x_not_y); // cannot be both true
      has_satisfy_x = has_satisfy_x_and_y 
                   || has_satisfy_x_not_y;
      if (!has_satisfy_x) {
        return true;
      }
      x_supp++;

      const bool kSatisfyXSuppCallbackRet
                = satisfy_x_supp_callback(pivoted_pattern_to_data_graph_match);
      if (has_satisfy_x_not_y) {
        return kSatisfyXSuppCallbackRet;
      }
      assert(has_satisfy_x_and_y);
      xy_supp++;
      const bool kSatisfyXYSuppCallbackRet
                = satisfy_xy_supp_callback(pivoted_pattern_to_data_graph_match);
      return kSatisfyXSuppCallbackRet
         && kSatisfyXYSuppCallbackRet;
    }
    /* ########################################################### *
     * ##  pivoted_pattern is the same as pivoted_rhs_pattern,  ## *
     * ##  directly whether the match satisfy y literals here   ## *
     * ########################################################### */
    assert(pivoted_rhs_pattern_to_data_graph_partial_match.size()
            == pivoted_pattern_to_data_graph_match.size());

    /* ########################################################### *
     * ##  get the paritial match from pivoted_literal_pattern  ## *
     * ##  to data_graph                                        ## *
     * ########################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
    pivoted_literal_pattern_to_pivoted_pattern_partial_match);
    if (!pivoted_rhs_gar_same_as_pivoted_literal_gar) {
      assert( GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_pattern));
      assert(!        SameGar    (pivoted_rhs_gar,     pivoted_literal_gar));
      assert(!GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_literal_pattern));
      
      assert(!has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph, 
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set,
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    }
    else {
      // pivoted_rhs_pattern is the same as pivoted_literal_pattern
      // directly match to pivoted_literal_pattern
      /* ############################################### *
       * ##  get the paritial match from gar_pattern  ## *
       * ##  to data_graph                            ## *
       * ############################################### */
      assert(        SameGar    (pivoted_rhs_gar, pivoted_literal_gar));
      assert(GUNDAM::SamePattern(pivoted_pattern, pivoted_literal_pattern));
      /* ######################################################## *
       * ##  chech whether the match satisfy lhs literal here  ## *
       * ######################################################## */
      bool satisify_x_literals = true;
      for (const auto& x_literal_ptr : pivoted_literal_gar.x_literal_set()) {
        assert(x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_partial_match));
        if (x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_partial_match)) {
          continue;
        }
        satisify_x_literals = false;
        break;
      }
      if (!satisify_x_literals) {
        // continue matching
        assert(!has_satisfy_x);
        return true;
      }
      if (!pivoted_literal_pattern_same_as_gar_pattern) {
        GUNDAM::Match<GraphPatternType, DataGraphType>
              gar_pattern_to_data_graph_partial_match
        = pivoted_pattern_to_data_graph_match(
              gar_pattern_to_pivoted_pattern_partial_match);
        assert(!has_satisfy_x);
        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                gar_pattern,   data_graph, 
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
      }
      else {
        assert(GUNDAM::SamePattern(pivoted_pattern, gar_pattern));
        // pivoted_pattern is the same as gar_pattern, no longer need to match
        assert(!has_satisfy_x);
        has_satisfy_x = true;
        auto data_graph_x_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
        auto data_graph_y_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

        support_candidate_set[data_graph_x_handle].insert(data_graph_y_handle);
      }
    }
    assert(!has_satisfy_x_and_y);
    assert(!has_satisfy_x_not_y);

    if (!has_satisfy_x) {
      return true;
    }
    x_supp++;

    GUNDAM::Match<GraphPatternType, DataGraphType>
          gar_pattern_to_data_graph_partial_match
    = pivoted_pattern_to_data_graph_match(
          gar_pattern_to_pivoted_pattern_partial_match);

    const bool kSatisfyXSuppCallbackRet
              = satisfy_x_supp_callback(gar_pattern_to_data_graph_partial_match);

    /* ######################################################## *
    * ##  chech whether the match satisfy rhs literal here  ## *
    * ######################################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_partial_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_partial_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }
    if (!satisify_y_literals) {
      return kSatisfyXSuppCallbackRet;
    }
    // satisfy y literals
    xy_supp++;
    // continue matching
    const bool kSatisfyXYSuppCallbackRet
              = satisfy_xy_supp_callback(gar_pattern_to_data_graph_partial_match);
    return kSatisfyXSuppCallbackRet
       && kSatisfyXYSuppCallbackRet;
  };

  MatchType pivoted_pattern_partial_match;
  for (auto pivoted_pattern_vertex_it = pivoted_pattern.VertexBegin();
           !pivoted_pattern_vertex_it.IsDone();
            pivoted_pattern_vertex_it++) {
    auto target_handle = gar_pattern_to_data_graph_partial_match.MapTo(pivoted_pattern_vertex_it->id());
    if (!target_handle) {
      continue;
    }
    pivoted_pattern_partial_match.AddMap(pivoted_pattern_vertex_it,
                                         target_handle);
  }

  // first match pivoted gar_pattern, then match rhs gar_pattern
  GUNDAM::MatchUsingMatch<match_semantics,
                          GUNDAM::MatchAlgorithm::kDagDp,
                          GUNDAM::MergeNecConfig::kNotMerge>(
                          pivoted_pattern, data_graph,
                          pivoted_pattern_partial_match,
                          pivoted_pattern_candidate_set,
                          pivoted_pattern_prune_callback,
                          pivoted_pattern_match_callback,
                          time_limit);
  assert(x_supp >= xy_supp);
  return std::pair(x_supp, xy_supp);
}



template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<std::vector<typename    DataGraphType::VertexCounterType>,
          std::vector<typename    DataGraphType::VertexCounterType>>
  GarSuppForSet(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& common_gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::vector<gar::LiteralInfo<GraphPatternType,
                         DataGraphType>>& new_literals_to_check,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
     std::function<bool(const unsigned, const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const unsigned, const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_xy_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  all_satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> all_satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {
  bool store_match = false;
  std::vector<std::vector<GUNDAM::Match<GraphPatternType, DataGraphType>>> match_vec_for_all;
  return GarSuppForSet(common_gar, data_graph, gar_pattern_to_data_graph_partial_match,
                      new_literals_to_check,
                      gar_pattern_candidate_set_for_vertexes_contained_in_literals,
                      gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                      store_match,
                      match_vec_for_all,
                      satisfy_x_supp_callback,
                      satisfy_xy_supp_callback,
                      all_satisfy_x_supp_callback,
                      all_satisfy_xy_supp_callback,
                      time_limit,
                      time_limit_per_supp);
}




/*  ################################################  *
 *  ##                                            ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##          pivoted_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##      pivoted_rhs_pattern_match_callback    ##  *
 *  ##      // currently this level is omitted    ##  *
 *  ##      // since it must be same with         ##  *
 *  ##      // pivoted_pattern now                ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##common_pivoted_literal_pattern_match_callback#  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##  pivoted_literal_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##              gar_pattern_match_callback    ##  *
 *  ##                                            ##  *
 *  ################################################  */

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<std::vector<typename    DataGraphType::VertexCounterType>,
          std::vector<typename    DataGraphType::VertexCounterType>>
  IncGarSuppForSet(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& common_gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::vector<gar::LiteralInfo<GraphPatternType,
                         DataGraphType>>& new_literals_to_check,
  std::map<typename GUNDAM::VertexHandle<DataGraphType>::type,
              std::map<typename GUNDAM::VertexHandle<DataGraphType>::type,
                        std::set<unsigned>>>& support_candidate_set,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
     std::function<bool(const unsigned, const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const unsigned, const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_xy_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  all_satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> all_satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;
  using VertexIDType = typename GUNDAM::VertexID<GraphPatternType>::type;

  assert(time_limit == -1.0
      || time_limit > 0);

  assert(time_limit_per_supp == -1.0
      || time_limit_per_supp > 0);

  assert(pivoted_vertex_set.empty() 
      || pivoted_vertex_set.size() > gar_pattern_to_data_graph_partial_match.size());

  if (// time_limit is set but time_limit_per_supp is not set
      (time_limit_per_supp == -1.0 && time_limit != -1.0)
      // both time_limit and time_limit_per_supp are set
      // but there are time_limit_per_supp > time_limit
   || (time_limit_per_supp != -1.0 && time_limit != -1.0
    && time_limit_per_supp > time_limit)) {
    time_limit_per_supp = time_limit;
  }
  // match_vec_for_all.resize(new_literals_to_check.size());
  auto& common_gar_pattern = common_gar.pattern();

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;


  /* ############################################################## *
   * ##   obtain common pattern contained in all literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
        common_pivoted_literal_gar = _gar_supp::LiteralGar(common_gar, pivoted_vertex_set);
  auto& common_pivoted_literal_pattern = common_pivoted_literal_gar.pattern();

  /* ############################################################## *
   * ##   obtain patterns contained in all literals and pivotes   ## *
   * ############################################################## */
  std::vector<GraphAssociationRule<GraphPatternType, 
                          DataGraphType>>
        gars_vec, pivoted_literal_gars_vec;
  // std::vector<GraphAssociationRule<GraphPatternType, 
  //                         DataGraphType>>
  //       pivoted_literal_gars_to_check_vec;  
  for (unsigned idx = 0; idx < new_literals_to_check.size(); idx++) {
    GraphAssociationRule<GraphPatternType, 
                        DataGraphType>
                    gar_to_check = common_gar;
    gar_to_check.AddX(new_literals_to_check[idx]);
    gars_vec.emplace_back(gar_to_check);
    GraphAssociationRule<GraphPatternType, 
                        DataGraphType>
          pivoted_literal_gar =
                  _gar_supp::LiteralGar(gar_to_check, pivoted_vertex_set);
    pivoted_literal_gars_vec.emplace_back(pivoted_literal_gar);
  }

  /* ############################################################## *
   * ##   obtain pattern contained in rhs literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
         pivoted_rhs_gar = _gar_supp::RhsGar(common_gar, pivoted_vertex_set);
  auto&  pivoted_rhs_pattern = pivoted_rhs_gar.pattern();
  //assert(pivoted_literal_pattern.CountVertex()
  //        >= pivoted_rhs_pattern.CountVertex());
  //assert(GUNDAM::SubGraphOf(pivoted_rhs_pattern,
  //                      pivoted_literal_pattern));

  /* ############################################# *
   * ##   obtain pattern contained in pivotes   ## *
   * ############################################# */
  GraphPatternType pivoted_pattern = pivoted_vertex_set.empty()?
                                     pivoted_rhs_pattern
                                   : GUNDAM::PreserveVertexSet(common_gar.pattern(),
                                     pivoted_vertex_set);
  assert(pivoted_rhs_pattern.CountVertex() 
          >= pivoted_pattern.CountVertex());
  assert(pivoted_vertex_set.empty()
      || pivoted_pattern.CountVertex() 
      == pivoted_vertex_set.size());
  assert(GUNDAM::SubGraphOf(pivoted_pattern,
                        pivoted_rhs_pattern));
  #ifndef NDEBUG
  for (auto map_it = gar_pattern_to_data_graph_partial_match.MapBegin();
           !map_it.IsDone();
            map_it++) {
    assert(pivoted_pattern.FindVertex(map_it->src_handle()->id()));
  }
  #endif // NDEBUG

  /* ########################################## *
   * ##  parameters to be used in callbacks  ## *
   * ########################################## */
  const bool pivoted_pattern_same_as_pivoted_rhs_pattern 
               = GUNDAM::SamePattern(pivoted_pattern,
                                     pivoted_rhs_pattern);
  if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
    std::cout << "pivoted_pattern should be same with pivoted_rhs_pattern" << std::endl;
    return std::pair(std::vector<typename    DataGraphType::VertexCounterType>(),
          std::vector<typename    DataGraphType::VertexCounterType>());
  }

  const bool pivoted_rhs_gar_same_as_common_pivoted_literal_gar
               = SameGar(pivoted_rhs_gar, common_pivoted_literal_gar);
  
  std::vector<unsigned> pivoted_literal_pattern_same_as_common_pivoted_literal_pattern;
  for (unsigned pivoted_literal_gar_idx = 0;
                pivoted_literal_gar_idx < pivoted_literal_gars_vec.size();
                pivoted_literal_gar_idx++) {
    if (GUNDAM::SamePattern(common_pivoted_literal_pattern,
                            pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern())) {
      pivoted_literal_pattern_same_as_common_pivoted_literal_pattern.push_back(
                    pivoted_literal_gar_idx);
    }
  }

  std::vector<bool> pivoted_literal_pattern_same_as_gar_pattern(
          new_literals_to_check.size(), false);
  for (unsigned idx = 0;
          idx < pivoted_literal_pattern_same_as_gar_pattern.size();
          idx++) {
    pivoted_literal_pattern_same_as_gar_pattern[idx]
            = GUNDAM::SamePattern(pivoted_literal_gars_vec[idx].pattern(),
                                  gars_vec[idx].pattern());
  }

  std::set<VertexIDType> common_vertex_id_set;
  for (auto vertex_it = common_pivoted_literal_pattern.VertexBegin();
          !vertex_it.IsDone();
          vertex_it++) {
    common_vertex_id_set.emplace(vertex_it->id());
  }

  std::map<VertexIDType, std::set<VertexIDType>> delta_vid_to_gar_idx;

  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    auto &pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
    for (auto vertex_it = pivoted_literal_pattern.VertexBegin();
              !vertex_it.IsDone();
              vertex_it++) {
      if (common_vertex_id_set.find(vertex_it->id()) == common_vertex_id_set.end()) {
        delta_vid_to_gar_idx[vertex_it->id()].insert(idx);
      }
    }
  }


  std::vector<GraphAssociationRule<GraphPatternType, 
                        DataGraphType>> delta_pivoted_literal_gars_vec;
  std::vector<std::pair<VertexIDType, unsigned>> delta_vid_to_size;
  std::vector<VertexIDType> delta_vid_vec;
  std::vector<bool> visited_gar_idx(pivoted_literal_gars_vec.size(), false);
  std::vector<std::set<unsigned>>
        delta_pivoted_literal_gar_to_pivoted_literal_gar_vec;
  std::vector<std::set<unsigned>>
        delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec;
  std::map<unsigned, unsigned> gar_to_corresponding_delta_pivoted_literal_gar;


  for (auto [vid, idx_set] : delta_vid_to_gar_idx) {
    delta_vid_to_size.emplace_back(vid, idx_set.size());
  }

  // bool delta_pivoted_literal_gar_empty = false;
  // if (delta_vid_to_size.empty()) {
  //   delta_pivoted_literal_gar_empty = true;
  // }

  while (!delta_vid_to_size.empty()) {
    std::sort(delta_vid_to_size.begin(), delta_vid_to_size.end(),
        [](const std::pair<VertexIDType, unsigned> &p1,
           const std::pair<VertexIDType, unsigned> &p2) {
             return p1.second > p2.second; });
    auto front_vid = delta_vid_to_size.front().first;
    //delta_vid_vec.emplace_back(front_id);

    //create delta_pivoted_literal_gar
    auto delta_vid_handle = common_gar_pattern.FindVertex(front_vid);
    if (!delta_vid_handle) {
      std::cout << "should have the handle" << std::endl;
    }
    auto delta_pivoted_vertex_set = pivoted_vertex_set;
    delta_pivoted_vertex_set.emplace_back(delta_vid_handle);
    GraphAssociationRule<GraphPatternType, 
                        DataGraphType>
          delta_pivoted_literal_gar =
                  _gar_supp::LiteralGar(common_gar, delta_pivoted_vertex_set);
    delta_pivoted_literal_gars_vec.emplace_back(delta_pivoted_literal_gar);
    auto &delta_pivoted_literal_pattern = delta_pivoted_literal_gar.pattern();

    std::set<unsigned> delta_pivoted_literal_gar_to_pivoted_literal_gar;
    std::set<unsigned> delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern;
    for (auto idx : delta_vid_to_gar_idx[front_vid]) {
      auto &pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
      if (visited_gar_idx[idx]) {
        continue;
//        std::cout << "a pivoted_literal_gar should not be visited twice" << std::endl;
//        return std::pair(std::vector<uint64_t>(), std::vector<uint64_t>());
      }
      visited_gar_idx[idx] = true;

      // update parental relationship of delta_pivoted_literal gar to pivoted_literal_gar
      if (GUNDAM::SamePattern(pivoted_literal_pattern,
                              delta_pivoted_literal_pattern)) {
        delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern.insert(idx);
      } else {
        delta_pivoted_literal_gar_to_pivoted_literal_gar.insert(idx);
      }

      gar_to_corresponding_delta_pivoted_literal_gar[idx]
              = delta_pivoted_literal_gar_to_pivoted_literal_gar_vec.size();

      // updated the sorted vector delta_vid_to_size
      for (auto vertex_it = pivoted_literal_pattern.VertexBegin();
                !vertex_it.IsDone();
                vertex_it++) {
        if (common_vertex_id_set.find(vertex_it->id()) == common_vertex_id_set.end()) {
          bool found_vid = false;
          for (auto vec_it = delta_vid_to_size.begin();
                    vec_it <  delta_vid_to_size.end();) {
            if (vec_it->first == vertex_it->id()) {
              found_vid = true;
              unsigned current_size = vec_it->second;
              if (current_size == 0) {
                std::cout << "current size should not be 0" << std::endl;
              }
              current_size--;
              if (current_size == 0) {
                vec_it = delta_vid_to_size.erase(vec_it);
                continue;
              }
              *vec_it = std::make_pair(vertex_it->id(), current_size);
              vec_it++;
            } else {
              vec_it++;
            }
            //std::cout << "here 2" << std::endl;
          }
          if (!found_vid) {
            std::cout << "should have found a vid in delta_vid_to_size" << std::endl;
            return std::pair(std::vector<uint64_t>(), std::vector<uint64_t>());
          }
        }
      }
    //std::cout << "here1" << std::endl;
    }
    delta_pivoted_literal_gar_to_pivoted_literal_gar_vec.
          emplace_back(delta_pivoted_literal_gar_to_pivoted_literal_gar);
    delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec.
          emplace_back(delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern);
  }

  unsigned to_check = 0;
  for (unsigned i = 0; i < delta_pivoted_literal_gar_to_pivoted_literal_gar_vec.size(); i++) {
    to_check += delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[i].size();
    to_check += delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[i].size();
  }
  to_check += pivoted_literal_pattern_same_as_common_pivoted_literal_pattern.size();

  if (to_check != new_literals_to_check.size()) {
    std::cout << "size of to check is not same with new literals to check" << std::endl;
  }
    //std::cout << "finish delta tree generation" << std::endl;
  

  /* ############################ *
   * ##  from pivoted_pattern  ## *
   * ############################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match(
    pivoted_rhs_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    common_pivoted_literal_pattern_to_pivoted_pattern_partial_match(
    common_pivoted_literal_pattern,   pivoted_pattern, "same_id_map");

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    delta_pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    delta_pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec.emplace_back(
      delta_pivoted_literal_gars_vec[idx].pattern(), pivoted_pattern, "same_id_map");
  }


  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec.emplace_back(
      pivoted_literal_gars_vec[idx].pattern(), pivoted_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_pivoted_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_pivoted_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, pivoted_pattern, "same_id_map");
  }
  


  /* ################################ *
   * ##  from pivoted_rhs_pattern  ## *
   * ################################ */

  GUNDAM::Match<GraphPatternType, GraphPatternType>
    common_pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec(
      common_pivoted_literal_pattern,   pivoted_rhs_pattern, "same_id_map");

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    delta_pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    delta_pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec.emplace_back(
      delta_pivoted_literal_gars_vec[idx].pattern(), pivoted_rhs_pattern, "same_id_map");
  }  

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match_vec.emplace_back(
      pivoted_literal_gars_vec[idx].pattern(), pivoted_rhs_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_pivoted_rhs_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_pivoted_rhs_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, pivoted_rhs_pattern, "same_id_map");
  }

  /* ########################################### *
   * ##  from common_pivoted_literal_pattern  ## *
   * ########################################### */
  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    delta_pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    delta_pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec.emplace_back(
      delta_pivoted_literal_gars_vec[idx].pattern(), common_pivoted_literal_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
    pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec.emplace_back(
      pivoted_literal_gars_vec[idx].pattern(), common_pivoted_literal_pattern, "same_id_map");
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, common_pivoted_literal_pattern, "same_id_map");
  }

  /* ########################################## *
   * ##  from delta_pivoted_literal_pattern  ## *
   * ########################################## */  
  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec(
            pivoted_literal_gars_vec.size());
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto pivoted_literal_gar_idx :
              delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[idx]) {
      pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                pivoted_literal_gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
                  pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern(),
                  delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }

  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto pivoted_literal_gar_idx :
              delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[idx]) {
      pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                pivoted_literal_gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
                  pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern(),
                  delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }

  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec(gars_vec.size());
  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto gar_idx : // gar_pattern_idx is the same with pivoted_literal_pattern_idx
              delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[idx]) {
      gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
            common_gar_pattern,
            delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }

  for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
    for (auto gar_idx :
              delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[idx]) {
      gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[gar_idx]
            = GUNDAM::Match<GraphPatternType, GraphPatternType>(
            common_gar_pattern,
            //gars_vec[gar_idx].pattern(),
            delta_pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
    }
  }


  /* #################################### *
   * ##  from pivoted_literal_pattern  ## *
   * #################################### */
  std::vector<GUNDAM::Match<GraphPatternType, GraphPatternType>>
    gar_pattern_to_pivoted_literal_pattern_partial_match_vec;
  for (unsigned idx = 0; idx < gars_vec.size(); idx++) {
    gar_pattern_to_pivoted_literal_pattern_partial_match_vec.emplace_back(
      common_gar_pattern, pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
      //gars_vec[idx].pattern(), pivoted_literal_gars_vec[idx].pattern(), "same_id_map");
  }

  const auto kRhsLiteralInfo
         = (*common_gar.y_literal_set()
                .begin())->info();

  auto rhs_literal_x_id = kRhsLiteralInfo.x_id();
  auto rhs_literal_y_id = kRhsLiteralInfo.y_id();

  // GUNDAM::Match<GraphPatternType, GraphPatternType>
  //   gar_pattern_to_pivoted_literal_pattern_partial_match(
  //   gar_pattern,   pivoted_literal_pattern, "same_id_map");


  //####################################
  //* prepare candiate set for new gars*
  //####################################
  CandidateSetContainer pivoted_pattern_candidate_set,
                    pivoted_rhs_pattern_candidate_set,
         common_pivoted_literal_pattern_candidate_set;

  std::vector<CandidateSetContainer> 
          delta_pivoted_literal_pattern_candidate_set(delta_pivoted_literal_gars_vec.size()),
          pivoted_literal_pattern_candidate_set(pivoted_literal_gars_vec.size());

  for (const auto& [gar_pattern_candidate_ptr,
                    gar_pattern_candidate]
                  : gar_pattern_candidate_set_for_vertexes_contained_in_literals) {
    assert(gar_pattern_candidate_ptr);
    if (!gar_pattern_candidate_ptr) {
        std::cout << "error candidate" << std::endl;
    return std::pair(std::vector<typename    DataGraphType::VertexCounterType>(),
          std::vector<typename    DataGraphType::VertexCounterType>());

    }
     auto kCommonPivotedLiteralPatternPtr 
              = common_pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
    if (kCommonPivotedLiteralPatternPtr) {
      // this vertex is contained in common_pivoted_literal_pattern
      auto [ common_pivoted_literal_pattern_candidate_set_it,
             common_pivoted_literal_pattern_candidate_set_ret ]
           = common_pivoted_literal_pattern_candidate_set.emplace(
            kCommonPivotedLiteralPatternPtr,
            gar_pattern_candidate);
      // should added successfully
      assert(common_pivoted_literal_pattern_candidate_set_ret);

      for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
        auto& delta_pivoted_literal_pattern = delta_pivoted_literal_gars_vec[idx].pattern();
        auto kDeltaPivotedLiteralPatternPtr
              = delta_pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        auto [delta_pivoted_literal_pattern_candidate_set_it,
              delta_pivoted_literal_pattern_candidate_set_ret ]
            = delta_pivoted_literal_pattern_candidate_set[idx].emplace(
              kDeltaPivotedLiteralPatternPtr,
              gar_pattern_candidate);
        assert(delta_pivoted_literal_pattern_candidate_set_ret);
        // should added successfully
      }

      for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
        auto& pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
        auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        auto [ pivoted_literal_pattern_candidate_set_it,
               pivoted_literal_pattern_candidate_set_ret ]
            = pivoted_literal_pattern_candidate_set[idx].emplace(
              kPivotedLiteralPatternPtr,
              gar_pattern_candidate);
        assert(pivoted_literal_pattern_candidate_set_ret);
        // should added successfully
      }

      const auto kPivotedRhsPatternPtr
                = pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id());
      if (kPivotedRhsPatternPtr) {
        // this vertex is contained in pivoted_rhs_pattern
        auto [ pivoted_rhs_pattern_candidate_set_it,
               pivoted_rhs_pattern_candidate_set_ret ]
             = pivoted_rhs_pattern_candidate_set.emplace(
              kPivotedRhsPatternPtr,
                gar_pattern_candidate);
        // should added successfully
        assert(pivoted_rhs_pattern_candidate_set_ret);
        const auto kPivotedPatternPtr
                  = pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedPatternPtr){
          // this vertex is contained in pivoted_rhs_pattern
          auto [ pivoted_pattern_candidate_set_it,
                 pivoted_pattern_candidate_set_ret ]
               = pivoted_pattern_candidate_set.emplace(
                 kPivotedPatternPtr,
                  gar_pattern_candidate);
          // should added successfully
          assert(pivoted_pattern_candidate_set_ret);
        }
      }
    } else {
      for (unsigned idx = 0; idx < delta_pivoted_literal_gars_vec.size(); idx++) {
        auto& delta_pivoted_literal_pattern = delta_pivoted_literal_gars_vec[idx].pattern();
        auto kDeltaPivotedLiteralPatternPtr
              = delta_pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kDeltaPivotedLiteralPatternPtr) {
          auto [ delta_pivoted_literal_pattern_candidate_set_it,
                delta_pivoted_literal_pattern_candidate_set_ret ]
              = delta_pivoted_literal_pattern_candidate_set[idx].emplace(
                kDeltaPivotedLiteralPatternPtr,
                gar_pattern_candidate);
          // should added successfully
          assert(delta_pivoted_literal_pattern_candidate_set_ret);
        }
      }

      for (unsigned idx = 0; idx < pivoted_literal_gars_vec.size(); idx++) {
        auto& pivoted_literal_pattern = pivoted_literal_gars_vec[idx].pattern();
        auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedLiteralPatternPtr) {
          auto [ pivoted_literal_pattern_candidate_set_it,
                pivoted_literal_pattern_candidate_set_ret ]
              = pivoted_literal_pattern_candidate_set[idx].emplace(
                kPivotedLiteralPatternPtr,
                gar_pattern_candidate);
          // should added successfully
          assert(pivoted_literal_pattern_candidate_set_ret);
        }
      }
    }
  }

  if (new_literals_to_check.size() == 0) {
    std::cout << "literal size is 0" << std::endl;
  }

  std::vector<bool> satisfy_x_callback_gar(gars_vec.size(), false),
                    satisfy_xy_callback_gar(gars_vec.size(), false),
                    satisfy_x_callback_delta_gar(
                        delta_pivoted_literal_gars_vec.size(), false);
  
  std::vector<typename DataGraphType::VertexCounterType>
                        x_supp_counter_vec(gars_vec.size(), 0),
                       xy_supp_counter_vec(gars_vec.size(), 0);

  std::vector<unsigned> satisfy_x_callback_delta_gar_num(
                              delta_pivoted_literal_gars_vec.size(), 0);

  std::vector<bool> has_x_match_gar(gars_vec.size(), false),
                    all_has_x_match_delta_pivoted_literal_gar(
                          delta_pivoted_literal_gars_vec.size(), false);
  
  // std::vector<bool> gar_satisfy_x_supp(gars_vec.size(), false),
  //                   gar_has_match_vec(gars_vec.size(), false);
  //                   delta_pivoted_literal_gar_has_match(
  //                         delta_pivoted_literal_gars_vec.size(), false),
                    
  unsigned gar_idx = 0,
           pivoted_literal_gar_idx = 0,
           delta_pivoted_literal_gar_idx = 0;
  
  bool all_has_satisfy_x = false;

  // bool has_satisfy_x       = false,
  //      has_satisfy_x_not_y = false,
  //      has_satisfy_x_and_y = false;

  // typename DataGraphType::VertexCounterType x_supp = 0,
  //                                          xy_supp = 0;

  /* ################################################ *
   * ##  callbacks from gar_pattern to data graph  ## *
   * ################################################ */
  std::function<bool(const MatchType&)> 
    gar_pattern_prune_callback
    = [](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // prune nothing, continue the matching
    return false;
  };

  // based on the partial match from literal_kPattern to data graph
  // complete the match of the entire gar_pattern 
  std::function<bool(const MatchType&)> 
    gar_pattern_match_callback 
    = [&all_has_satisfy_x,
       &has_x_match_gar,
       &gar_idx,
       &x_supp_counter_vec,
       &rhs_literal_x_id,
       &rhs_literal_y_id,
       &support_candidate_set](const MatchType& 
    gar_pattern_to_data_graph_match) -> bool {
    // only need to complete the match, find one match that 
    // satisfy all literal is enough
    all_has_satisfy_x = true;
    x_supp_counter_vec[gar_idx]++;
    has_x_match_gar[gar_idx] = true;

    // if (store_match) {
    //   match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_match);
    // }
    auto data_graph_x_handle = gar_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
    auto data_graph_y_handle = gar_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

    support_candidate_set[data_graph_x_handle][data_graph_y_handle].insert(gar_idx);    

    // terminate matching
    return false;
  };

  /* ############################################################ *
   * ##  callbacks from pivoted_literal_pattern to data graph  ## *
   * ############################################################ */
  // if the current pivoted_literal_pattern_to_data_graph_match
  // has already violate a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr :
            pivoted_literal_gars_vec[pivoted_literal_gar_idx].x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_literal_pattern_to_data_graph_match)){
        // this literal is not considered in the current
        // pivoted_literal_pattern_to_data_graph_match,
        // this match cannot be pruned
        //
        // move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_literal_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        //std::cout << "exist x not match pivoted" << std::endl;
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };


  std::function<bool(const MatchType&)>
    pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    pivoted_literal_pattern_to_data_graph_match) -> bool {
    all_has_satisfy_x = false;
    gar_idx = pivoted_literal_gar_idx;
    assert(!has_x_match_gar[gar_idx]);

    if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
      // pivoted_rhs_gar is not the same as pivoted_literal_gar
      // first match the pivoted_rhs_literal, then match
      // the pivoted_literal_gar

      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
        = pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_pivoted_literal_pattern_partial_match_vec[
                            gar_idx]);
      assert(!all_has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              common_gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
      if (all_has_satisfy_x) {
        assert(has_x_match_gar[gar_idx]);
        //has_x_match_gar[gar_idx] = true;
        // has found satisfy X, terminate matching
        return false;
      }
      // continue matching
      return true;
    }
    assert(!all_has_satisfy_x);
    // gar_pattern of literal_gar is the same as the
    // entire gar_pattern, does not need further match
    all_has_satisfy_x = true;

    // if (store_match) {
    //   GUNDAM::Match<GraphPatternType, DataGraphType>
    //                   gar_pattern_to_data_graph_entire_match 
    //     = pivoted_literal_pattern_to_data_graph_match(
    //                   gar_pattern_to_pivoted_literal_pattern_partial_match_vec[
    //                         gar_idx]);
    //   match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
    // }

    auto data_graph_x_handle = pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
    auto data_graph_y_handle = pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

    support_candidate_set[data_graph_x_handle][data_graph_y_handle].insert(gar_idx); 


    x_supp_counter_vec[gar_idx]++;
    has_x_match_gar[gar_idx] = true;
    // has found satisfy X, terminate matching
    return false;
  };

  std::function<bool(const MatchType&)>
    delta_pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    delta_pivoted_literal_pattern_to_data_graph_match) -> bool {
    //prune nothing
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    delta_pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    delta_pivoted_literal_pattern_to_data_graph_match) -> bool {
    unsigned gar_num_no_need_check = 0;

    if (all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]) {
      std::cout << "should not be here delta" << std::endl;
    }

    assert(!all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]);

    for (auto current_pivoted_literal_gar_idx :
            delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[
                  delta_pivoted_literal_gar_idx]) {
      pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
      gar_idx = pivoted_literal_gar_idx;

      if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
        gar_num_no_need_check++;
        continue;
      }
      all_has_satisfy_x = false;
      
      if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
                  delta_pivoted_literal_pattern_to_data_graph_match)) {
        continue;
      }

      if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
        GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
          = delta_pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                          gar_idx]);

        GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              common_gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
        if (all_has_satisfy_x) {
          gar_num_no_need_check++;
        }
      } else {
        // continue matching
        // if (store_match) {
        //   GUNDAM::Match<GraphPatternType, DataGraphType>
        //                   gar_pattern_to_data_graph_entire_match 
        //     = delta_pivoted_literal_pattern_to_data_graph_match(
        //               gar_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
        //                   gar_idx]);
        //   match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
        // }

        auto data_graph_x_handle = delta_pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
        auto data_graph_y_handle = delta_pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

        support_candidate_set[data_graph_x_handle][data_graph_y_handle].insert(gar_idx); 
        x_supp_counter_vec[gar_idx]++;
        has_x_match_gar[gar_idx] = true;
        gar_num_no_need_check++;
      }
    }

    for (auto current_pivoted_literal_gar_idx :
              delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[
                  delta_pivoted_literal_gar_idx]) {
      pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
      gar_idx = pivoted_literal_gar_idx;
      if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
        gar_num_no_need_check++;
        continue;
      }

      all_has_satisfy_x = false;
      auto& pivoted_literal_pattern
              = pivoted_literal_gars_vec[pivoted_literal_gar_idx].pattern();
      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      pivoted_literal_pattern_to_data_graph_partial_match
        = delta_pivoted_literal_pattern_to_data_graph_match(
                pivoted_literal_pattern_to_delta_pivoted_literal_pattern_partial_match_vec[
                    pivoted_literal_gar_idx]);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_literal_pattern,   data_graph,
                              pivoted_literal_pattern_to_data_graph_partial_match,
                              pivoted_literal_pattern_candidate_set[pivoted_literal_gar_idx],
                              pivoted_literal_pattern_prune_callback,
                              pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
      if (all_has_satisfy_x) {
        gar_num_no_need_check++;
      }
    }

    all_has_satisfy_x = false;
    if (gar_num_no_need_check
                  == (delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[
                      delta_pivoted_literal_gar_idx].size()
                  + delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[
                      delta_pivoted_literal_gar_idx].size())) {
      all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx] = true;
      all_has_satisfy_x = true;
      return false;
    }
    return true;
  };


  /* ############################################################ *
   * ##  callbacks from common_pivoted_literal_pattern to data graph  ## *
   * ############################################################ */
  // if the current pivoted_literal_pattern_to_data_graph_match
  // has already violate a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    common_pivoted_literal_pattern_prune_callback
    = [&](const MatchType& 
    common_pivoted_literal_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : common_pivoted_literal_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(common_pivoted_literal_pattern_to_data_graph_match)){
        // this literal is not considered in the current
        // pivoted_literal_pattern_to_data_graph_match,
        // this match cannot be pruned
        //
        // move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(common_pivoted_literal_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        //std::cout << "exist x not match common" << std::endl;
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };


  std::function<bool(const MatchType&)>
    common_pivoted_literal_pattern_match_callback
    = [&](const MatchType& 
    common_pivoted_literal_pattern_to_data_graph_match) -> bool {
    unsigned gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback = 0;
    for (auto current_pivoted_literal_gar_idx :
            pivoted_literal_pattern_same_as_common_pivoted_literal_pattern) {
      pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
      gar_idx = pivoted_literal_gar_idx;

      if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
        continue;
      }

      all_has_satisfy_x = false;

      if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
            common_pivoted_literal_pattern_to_data_graph_match)) {
        continue;
      }

      if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
        GUNDAM::Match<GraphPatternType, DataGraphType>  
                      gar_pattern_to_data_graph_partial_match
          = common_pivoted_literal_pattern_to_data_graph_match(
                      gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec[
                          gar_idx]);

        GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kMerge>(
                              common_gar_pattern,   data_graph,
                              gar_pattern_to_data_graph_partial_match,
                              gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                              gar_pattern_prune_callback,
                              gar_pattern_match_callback,
                              time_limit_per_supp);
        if (all_has_satisfy_x) {
          gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
        }
      } else {
        // if (store_match) {
        //   GUNDAM::Match<GraphPatternType, DataGraphType>
        //                   gar_pattern_to_data_graph_entire_match 
        //     = common_pivoted_literal_pattern_to_data_graph_match(
        //               gar_pattern_to_common_pivoted_literal_pattern_partial_match_vec[
        //                   gar_idx]);
        //   match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
        // }

        auto data_graph_x_handle = common_pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
        auto data_graph_y_handle = common_pivoted_literal_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

        support_candidate_set[data_graph_x_handle][data_graph_y_handle].insert(gar_idx); 

        has_x_match_gar[gar_idx] = true;
        x_supp_counter_vec[gar_idx]++;
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
      }
    }

    for (delta_pivoted_literal_gar_idx = 0;
         delta_pivoted_literal_gar_idx < delta_pivoted_literal_gars_vec.size();
         delta_pivoted_literal_gar_idx++) {
    // delta_pivoted_literal_gar can not be the same with common_pivoted_literal_gar

      if (satisfy_x_callback_delta_gar[delta_pivoted_literal_gar_idx]
              || all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]) {
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
        continue;
      }
      auto& delta_pivoted_literal_pattern
              = delta_pivoted_literal_gars_vec[delta_pivoted_literal_gar_idx].pattern();
      GUNDAM::Match<GraphPatternType, DataGraphType>  
                      delta_pivoted_literal_pattern_to_data_graph_partial_match
        = common_pivoted_literal_pattern_to_data_graph_match(
                delta_pivoted_literal_pattern_to_common_pivoted_literal_pattern_partial_match_vec[
                    delta_pivoted_literal_gar_idx]);
      all_has_satisfy_x = false;
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              delta_pivoted_literal_pattern,   data_graph,
                              delta_pivoted_literal_pattern_to_data_graph_partial_match,
                              delta_pivoted_literal_pattern_candidate_set[delta_pivoted_literal_gar_idx],
                              delta_pivoted_literal_pattern_prune_callback,
                              delta_pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
      if (all_has_satisfy_x) {
        gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback++;
      }
    }
    
    all_has_satisfy_x = false;
    if (gar_num_no_need_check_in_common_pivoted_literal_pattern_match_callback
            == delta_pivoted_literal_gars_vec.size()
              + pivoted_literal_pattern_same_as_common_pivoted_literal_pattern.size()) {
      all_has_satisfy_x = true;
      return false;
    }
    return true;
  };



  /* ######################################################## *
   * ##  callbacks from pivoted_rhs_pattern to data graph  ## *
   * ######################################################## */
  // if the current match_state has already violate
  // a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_prune_callback
    = [&](const MatchType&
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    std::cout << "should not use this prune callack" << std::endl;
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_match_callback
    = [&](const MatchType& 
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    std::cout << "should not use this match callback" << std::endl;
    return true;
  };

  /* #################################################### *
   * ##  callbacks from pivoted_pattern to data graph  ## *
   * #################################################### */
  std::function<bool(const MatchType&)>
    pivoted_pattern_prune_callback
    = [&](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    // prune nothing
    auto data_graph_x_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
    
    auto data_graph_y_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);

    if ((!data_graph_y_handle) || (!data_graph_x_handle)) {
      return false;
    }
    if ((support_candidate_set.find(data_graph_x_handle) != support_candidate_set.end())
        && (support_candidate_set[data_graph_x_handle].find(data_graph_y_handle)
              != support_candidate_set[data_graph_x_handle].end())) {
      // the pivot pair has already been verified as legal
      if (support_candidate_set[data_graph_x_handle][data_graph_y_handle].size()
            == gars_vec.size()) {
        return true;
      }
    }

    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_pattern_match_callback
    = [&](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    /* ######################################## *
     * ##  get a match for pivot vertex set  ## *
     * ######################################## */
    all_has_satisfy_x       = false;
    has_x_match_gar.clear();
    has_x_match_gar.resize(gars_vec.size(), false);
    all_has_x_match_delta_pivoted_literal_gar.clear();
    all_has_x_match_delta_pivoted_literal_gar.resize(
              delta_pivoted_literal_gars_vec.size(), false);

    {
      auto data_graph_x_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
      
      auto data_graph_y_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);

      if ((!data_graph_y_handle) || (!data_graph_x_handle)) {
        return false;
      }
      if ((support_candidate_set.find(data_graph_x_handle) != support_candidate_set.end())
          && (support_candidate_set[data_graph_x_handle].find(data_graph_y_handle)
                != support_candidate_set[data_graph_x_handle].end())) {
        // the pivot pair has already been verified as legal
        for (auto id : support_candidate_set[data_graph_x_handle][data_graph_y_handle]) {
          has_x_match_gar[id] = true;
        }
      }
    }


    /* ####################################################### *
     * ##  get the paritial match from pivoted_rhs_pattern  ## *
     * ##  to data_graph                                    ## *
     * ####################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_rhs_pattern_to_data_graph_partial_match
        = pivoted_pattern_to_data_graph_match(
          pivoted_rhs_pattern_to_pivoted_pattern_partial_match);

    if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar
      std::cout << "it has not been implemented and should not be here" << std::endl;
      return true;
    }
    /* ########################################################### *
     * ##  pivoted_pattern is the same as pivoted_rhs_pattern,  ## *
     * ##  directly whether the match satisfy y literals here   ## *
     * ########################################################### */
    assert(pivoted_rhs_pattern_to_data_graph_partial_match.size()
            == pivoted_pattern_to_data_graph_match.size());

    /* ########################################################### *
     * ##  get the paritial match from pivoted_literal_pattern  ## *
     * ##  to data_graph                                        ## *
     * ########################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    common_pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
            common_pivoted_literal_pattern_to_pivoted_pattern_partial_match);


    if (!pivoted_rhs_gar_same_as_common_pivoted_literal_gar) {
      //std::cout << "rhs not same as common" << std::endl;
      assert( GUNDAM::SamePattern(pivoted_rhs_pattern, pivoted_pattern));
      assert(!        SameGar    (pivoted_rhs_gar,     common_pivoted_literal_gar));
      assert(!GUNDAM::SamePattern(pivoted_rhs_pattern, common_pivoted_literal_pattern));

      assert(!all_has_satisfy_x);
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              common_pivoted_literal_pattern,   data_graph, 
                              common_pivoted_literal_pattern_to_data_graph_partial_match,
                              common_pivoted_literal_pattern_candidate_set,
                              common_pivoted_literal_pattern_prune_callback,
                              common_pivoted_literal_pattern_match_callback,
                              time_limit_per_supp);
    } else {
      //pivoted_rhs_gar is the same as common_pivoted_literal_gar
      assert(        SameGar    (pivoted_rhs_gar, common_pivoted_literal_gar));
      assert(GUNDAM::SamePattern(pivoted_pattern, common_pivoted_literal_pattern));

      // verify x_literals for common_pivoted_literal_gar
      bool satisify_x_literals = true;
      for (const auto& x_literal_ptr : common_pivoted_literal_gar.x_literal_set()) {
        assert(x_literal_ptr->MappedBy(common_pivoted_literal_pattern_to_data_graph_partial_match));
        if (x_literal_ptr->Satisfy(common_pivoted_literal_pattern_to_data_graph_partial_match)) {
          continue;
        }
        satisify_x_literals = false;
        break;
      }
      if (!satisify_x_literals) {
        // continue matching
        //std::cout << "exist x not match" << std::endl;
        assert(!satisify_x_literals);
        return true;
      }


      for (auto current_pivoted_literal_gar_idx :
              pivoted_literal_pattern_same_as_common_pivoted_literal_pattern) {
        pivoted_literal_gar_idx = current_pivoted_literal_gar_idx;
        gar_idx = pivoted_literal_gar_idx;

        if (satisfy_x_callback_gar[gar_idx] || has_x_match_gar[gar_idx]) {
          continue;
        }

        // GUNDAM::Match<GraphPatternType, DataGraphType>  
        //                 pivoted_literal_pattern_to_data_graph_partial_match
        //   = pivoted_pattern_to_data_graph_match(
        //             pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec[
        //                 pivoted_literal_gar_idx]);

        // if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
        //           pivoted_literal_pattern_to_data_graph_partial_match)) {
        //   continue;
        // }


        if (!new_literals_to_check[pivoted_literal_gar_idx].Satisfy(
                  pivoted_pattern_to_data_graph_match)) {
          continue;
        }

        if (!pivoted_literal_pattern_same_as_gar_pattern[pivoted_literal_gar_idx]) {
          GUNDAM::Match<GraphPatternType, DataGraphType>  
                        gar_pattern_to_data_graph_partial_match
            = pivoted_pattern_to_data_graph_match(
                        gar_pattern_to_pivoted_pattern_partial_match_vec[
                            gar_idx]);

          GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kMerge>(
                                common_gar_pattern,   data_graph,
                                gar_pattern_to_data_graph_partial_match,
                                gar_pattern_candidate_set_for_vertexes_not_contained_in_literals,
                                gar_pattern_prune_callback,
                                gar_pattern_match_callback,
                                time_limit_per_supp);
        } else {
          auto data_graph_x_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_x_id);
          auto data_graph_y_handle = pivoted_pattern_to_data_graph_match.MapTo(rhs_literal_y_id);    

          support_candidate_set[data_graph_x_handle][data_graph_y_handle].insert(gar_idx);           
          // if (store_match) {
          //   GUNDAM::Match<GraphPatternType, DataGraphType>
          //                   gar_pattern_to_data_graph_entire_match 
          //     = pivoted_pattern_to_data_graph_match(
          //               gar_pattern_to_pivoted_pattern_partial_match_vec[
          //                   gar_idx]);
          //   match_vec_for_all[gar_idx].emplace_back(gar_pattern_to_data_graph_entire_match);
          // }

          has_x_match_gar[gar_idx] = true;
          x_supp_counter_vec[gar_idx]++;
        }
      }


      for (delta_pivoted_literal_gar_idx = 0;
          delta_pivoted_literal_gar_idx < delta_pivoted_literal_gars_vec.size();
          delta_pivoted_literal_gar_idx++) {
      // delta_pivoted_literal_gar can not be the same with common_pivoted_literal_gar

        if (satisfy_x_callback_delta_gar[delta_pivoted_literal_gar_idx]
              || all_has_x_match_delta_pivoted_literal_gar[delta_pivoted_literal_gar_idx]) {
          continue;
        }
        auto& delta_pivoted_literal_pattern
                = delta_pivoted_literal_gars_vec[delta_pivoted_literal_gar_idx].pattern();
        GUNDAM::Match<GraphPatternType, DataGraphType>  
                        delta_pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
                  delta_pivoted_literal_pattern_to_pivoted_pattern_partial_match_vec[
                      delta_pivoted_literal_gar_idx]);
        GUNDAM::MatchUsingMatch<match_semantics,
                                GUNDAM::MatchAlgorithm::kDagDp,
                                GUNDAM::MergeNecConfig::kNotMerge>(
                                delta_pivoted_literal_pattern,   data_graph,
                                delta_pivoted_literal_pattern_to_data_graph_partial_match,
                                delta_pivoted_literal_pattern_candidate_set[delta_pivoted_literal_gar_idx],
                                delta_pivoted_literal_pattern_prune_callback,
                                delta_pivoted_literal_pattern_match_callback,
                                time_limit_per_supp);
      }
    }

     GUNDAM::Match<GraphPatternType, DataGraphType>
           gar_pattern_to_data_graph_partial_match
     = pivoted_pattern_to_data_graph_match(
           gar_pattern_to_pivoted_pattern_partial_match_vec[0]);


    bool exist_x_satisfy = false;
    for (unsigned gars_vec_idx = 0;
                  gars_vec_idx < gars_vec.size();
                  gars_vec_idx++) {
      if (has_x_match_gar[gars_vec_idx]) {
        exist_x_satisfy = true;

        const bool kSatisfyXSuppCallbackRet
              = satisfy_x_supp_callback(gars_vec_idx, gar_pattern_to_data_graph_partial_match);
        if (!kSatisfyXSuppCallbackRet) {
          satisfy_x_callback_gar[gars_vec_idx] = true;
        }
      }
    }

    if (!exist_x_satisfy) {
      return true;
    }

    /* ######################################################## *
    * ##  chech whether the match satisfy rhs literal here  ## *
    * ######################################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_partial_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_partial_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }
    const bool kAllSatisfyXSuppCallbackRet
          = all_satisfy_x_supp_callback(gar_pattern_to_data_graph_partial_match);
    if (!satisify_y_literals) {
      return kAllSatisfyXSuppCallbackRet;
    }
    // satisfy y literals

    for (unsigned gars_vec_idx = 0;
                  gars_vec_idx < gars_vec.size();
                  gars_vec_idx++) {
      if (has_x_match_gar[gars_vec_idx]) {
        xy_supp_counter_vec[gars_vec_idx]++;
        const bool kSatisfyXYSuppCallbackRet
              = satisfy_xy_supp_callback(gars_vec_idx, gar_pattern_to_data_graph_partial_match);
        if (!kSatisfyXYSuppCallbackRet) {
          satisfy_x_callback_gar[gars_vec_idx] = true;
          if (gar_to_corresponding_delta_pivoted_literal_gar.find(gars_vec_idx)
                  == gar_to_corresponding_delta_pivoted_literal_gar.end()) {
            continue;
          }
          unsigned corresponding_delta_pivoted_literal_gar_idx
                = gar_to_corresponding_delta_pivoted_literal_gar[gars_vec_idx];
          satisfy_x_callback_delta_gar_num[corresponding_delta_pivoted_literal_gar_idx]++;
          if (satisfy_x_callback_delta_gar_num[corresponding_delta_pivoted_literal_gar_idx]
                < delta_pivoted_literal_gar_to_pivoted_literal_gar_vec[
                      corresponding_delta_pivoted_literal_gar_idx].size()
                  + delta_pivoted_literal_pattern_same_as_pivoted_literal_pattern_vec[
                      corresponding_delta_pivoted_literal_gar_idx].size()) {
            continue;
          }
          satisfy_x_callback_delta_gar[corresponding_delta_pivoted_literal_gar_idx] = true;
        }
      }
    }
    // continue matching
    const bool kAllSatisfyXYSuppCallbackRet
              = all_satisfy_xy_supp_callback(gar_pattern_to_data_graph_partial_match);
    return kAllSatisfyXYSuppCallbackRet
            && kAllSatisfyXSuppCallbackRet;

  };

  MatchType pivoted_pattern_partial_match;
  for (auto pivoted_pattern_vertex_it = pivoted_pattern.VertexBegin();
           !pivoted_pattern_vertex_it.IsDone();
            pivoted_pattern_vertex_it++) {
    auto target_handle = gar_pattern_to_data_graph_partial_match.MapTo(pivoted_pattern_vertex_it->id());
    if (!target_handle) {
      continue;
    }
    pivoted_pattern_partial_match.AddMap(pivoted_pattern_vertex_it,
                                         target_handle);
  }

  // first match pivoted gar_pattern, then match rhs gar_pattern
  GUNDAM::MatchUsingMatch<match_semantics,
                          GUNDAM::MatchAlgorithm::kDagDp,
                          GUNDAM::MergeNecConfig::kNotMerge>(
                          pivoted_pattern, data_graph,
                          pivoted_pattern_partial_match,
                          pivoted_pattern_candidate_set,
                          pivoted_pattern_prune_callback,
                          pivoted_pattern_match_callback,
                          time_limit);
  //assert(x_supp >= xy_supp);
  return std::pair(x_supp_counter_vec, xy_supp_counter_vec);
}




template<typename GraphPatternType,
         typename DataGraphType>
bool TransferMatchToID(GraphAssociationRule<GraphPatternType,
                                            DataGraphType>& gar,
                                            DataGraphType& data_graph,
                  std::vector<GUNDAM::Match<GraphPatternType,
                                            DataGraphType>>& match_vec,
    std::vector<std::pair<std::vector<typename GUNDAM::VertexID<DataGraphType>::type>,
                          std::vector<typename GUNDAM::EdgeID<DataGraphType>::type>>>&
                                            match_id_vec) {
  auto &pattern = gar.pattern();

  const auto y_literal_ptr = *(gar.y_literal_set().begin());
  if (gar.y_literal_set().Count() != 1) {
    std::cout << "y literal set size is wrong" << std::endl;
    return false;
  }
  const auto kRhsLiteralInfo
         = y_literal_ptr->info();
  auto rhs_literal_x_id = kRhsLiteralInfo.x_id();
  auto rhs_literal_y_id = kRhsLiteralInfo.y_id();

  for (auto match : match_vec) {
    auto data_graph_x_id = (match.MapTo(rhs_literal_x_id))->id();
    auto data_graph_y_id = (match.MapTo(rhs_literal_y_id))->id();
    std::vector<typename GUNDAM::VertexID<DataGraphType>::type> pivot_vertex_id;
    pivot_vertex_id.emplace_back(data_graph_x_id);
    pivot_vertex_id.emplace_back(data_graph_y_id);

    std::vector<typename GUNDAM::EdgeID<DataGraphType>::type> edge_id;
    for (auto pattern_vertex_it = pattern.VertexBegin();
             !pattern_vertex_it.IsDone();
              pattern_vertex_it++) {
      auto pattern_vertex_id = pattern_vertex_it->id();
      auto data_vertex_handle = match.MapTo(pattern_vertex_id);
      if (!data_vertex_handle) {
        std::cout << "error handle " << std::endl;
      }
      auto data_vertex_id = data_vertex_handle->id();

      for (auto out_edge_it = pattern_vertex_it->OutEdgeBegin();
               !out_edge_it.IsDone();
                out_edge_it++) {
        auto pattern_dst_id = out_edge_it->dst_handle()->id();
        auto data_dst_handle = match.MapTo(pattern_dst_id);
        if (!data_dst_handle) {
          std::cout << "error handle " << std::endl;
        }
        auto data_dst_id = data_dst_handle->id();

        auto out_edge_label = out_edge_it->label();
        auto data_out_edge_it = data_vertex_handle->OutEdgeBegin(out_edge_label, data_dst_handle);
        if (!data_out_edge_it) {
          std::cout << "not exist the certain edge" << std::endl;
          return false;
        }
        auto data_out_edge_id = data_out_edge_it->id();
        edge_id.emplace_back(data_out_edge_id);
      }
    }
    match_id_vec.emplace_back(std::pair(pivot_vertex_id, edge_id));
  }
  return true;
}

/*  ################################################  *
 *  ##                                            ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##          pivoted_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##      pivoted_rhs_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##  pivoted_literal_pattern_match_callback    ##  *
 *  ##                            |               ##  *
 *  ##                            |               ##  *
 *  ##                            V               ##  *
 *  ##              gar_pattern_match_callback    ##  *
 *  ##                                            ##  *
 *  ################################################  */

// return std::pair<x_supp, xy_supp>
template <enum GUNDAM::MatchSemantics match_semantics 
             = GUNDAM::MatchSemantics::kIsomorphism,
          typename GraphPatternType,
          typename    DataGraphType>
std::pair<typename    DataGraphType::VertexCounterType,
          typename    DataGraphType::VertexCounterType> 
  RhsGarSupp(GraphAssociationRule<GraphPatternType,
                                  DataGraphType>& gar,
             DataGraphType& data_graph,
  const GUNDAM::Match<GraphPatternType,
                         DataGraphType>& gar_pattern_to_data_graph_partial_match,
  const std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
     std::vector<typename GUNDAM::VertexHandle<   DataGraphType>::type>>& gar_pattern_candidate_set_for_vertexes_contained_in_literals,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)>  satisfy_x_supp_callback,
     std::function<bool(const GUNDAM::Match<GraphPatternType,
                                               DataGraphType>&)> satisfy_xy_supp_callback,
            double time_limit = -1.0,
            double time_limit_per_supp = -1.0,
  const std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>& pivoted_vertex_set 
      = std::vector<typename GUNDAM::VertexHandle<GraphPatternType>::type>()) {

  using GraphPatternVertexHandle = typename GUNDAM::VertexHandle<GraphPatternType>::type;
  using    DataGraphVertexHandle = typename GUNDAM::VertexHandle<   DataGraphType>::type;

  assert(time_limit == -1.0
      || time_limit > 0);

  assert(time_limit_per_supp == -1.0
      || time_limit_per_supp > 0);

  assert(pivoted_vertex_set.empty() 
      || pivoted_vertex_set.size() > gar_pattern_to_data_graph_partial_match.size());

  if (// time_limit is set but time_limit_per_supp is not set
      (time_limit_per_supp == -1.0 && time_limit != -1.0)
      // both time_limit and time_limit_per_supp are set
      // but there are time_limit_per_supp > time_limit
   || (time_limit_per_supp != -1.0 && time_limit != -1.0
    && time_limit_per_supp > time_limit)) {
    time_limit_per_supp = time_limit;
  }

  auto& gar_pattern = gar.pattern();

  using MatchType = GUNDAM::Match<GraphPatternType, 
                                     DataGraphType>;

  using CandidateSetContainer = std::map<GraphPatternVertexHandle, 
                                std::vector<DataGraphVertexHandle>>;

  /* ############################################################## *
   * ##   obtain pattern contained in all literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
        pivoted_literal_gar = _gar_supp::LiteralGar(gar, pivoted_vertex_set);
  auto& pivoted_literal_pattern = pivoted_literal_gar.pattern();

  /* ############################################################## *
   * ##   obtain pattern contained in rhs literals and pivotes   ## *
   * ############################################################## */
  GraphAssociationRule<GraphPatternType, 
                          DataGraphType> 
         pivoted_rhs_gar = _gar_supp::RhsGar(gar, pivoted_vertex_set);
  auto&  pivoted_rhs_pattern = pivoted_rhs_gar.pattern();
  assert(pivoted_literal_pattern.CountVertex()
          >= pivoted_rhs_pattern.CountVertex());
  assert(GUNDAM::SubGraphOf(pivoted_rhs_pattern,
                        pivoted_literal_pattern));

  /* ############################################# *
   * ##   obtain pattern contained in pivotes   ## *
   * ############################################# */
  GraphPatternType pivoted_pattern = pivoted_vertex_set.empty()?
                                     pivoted_rhs_pattern
                                   : GUNDAM::PreserveVertexSet(gar.pattern(),
                                     pivoted_vertex_set);
  assert(pivoted_rhs_pattern.CountVertex() 
          >= pivoted_pattern.CountVertex());
  assert(pivoted_vertex_set.empty()
      || pivoted_pattern.CountVertex() 
      == pivoted_vertex_set.size());
  assert(GUNDAM::SubGraphOf(pivoted_pattern,
                        pivoted_rhs_pattern));
  #ifndef NDEBUG
  for (auto map_it = gar_pattern_to_data_graph_partial_match.MapBegin();
           !map_it.IsDone();
            map_it++) {
    assert(pivoted_pattern.FindVertex(map_it->src_handle()->id()));
  }
  #endif // NDEBUG

  /* ########################################## *
   * ##  parameters to be used in callbacks  ## *
   * ########################################## */
  const bool pivoted_pattern_same_as_pivoted_rhs_pattern 
               = GUNDAM::SamePattern(pivoted_pattern,
                                     pivoted_rhs_pattern);

  const bool pivoted_rhs_gar_same_as_pivoted_literal_gar 
                           = SameGar(pivoted_rhs_gar,
                                     pivoted_literal_gar);

  const bool pivoted_literal_pattern_same_as_gar_pattern
           = GUNDAM::SamePattern(pivoted_literal_pattern,
                                             gar_pattern);

  /* ############################ *
   * ##  from pivoted_pattern  ## *
   * ############################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match(
    pivoted_rhs_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_pattern_partial_match(
    gar_pattern,   pivoted_pattern, "same_id_map");

  /* ################################ *
   * ##  from pivoted_rhs_pattern  ## *
   * ################################ */
  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    pivoted_literal_pattern_to_pivoted_rhs_pattern_partial_match(
    pivoted_literal_pattern,   pivoted_rhs_pattern, "same_id_map");

  GUNDAM::Match<GraphPatternType, GraphPatternType> 
    gar_pattern_to_pivoted_rhs_pattern_partial_match(
    gar_pattern,   pivoted_rhs_pattern, "same_id_map");

  /* #################################### *
   * ##  from pivoted_literal_pattern  ## *
   * #################################### */
  GUNDAM::Match<GraphPatternType, GraphPatternType>
    gar_pattern_to_pivoted_literal_pattern_partial_match(
    gar_pattern,   pivoted_literal_pattern, "same_id_map");

  CandidateSetContainer pivoted_pattern_candidate_set,
                    pivoted_rhs_pattern_candidate_set,
                pivoted_literal_pattern_candidate_set;

  // process pivoted_literal_pattern_candidate_set
  //         and pivoted_rhs_pattern_candidate_set
  for (const auto& [gar_pattern_candidate_ptr,
                    gar_pattern_candidate]
                  : gar_pattern_candidate_set_for_vertexes_contained_in_literals) {
    // should not be null
    assert(gar_pattern_candidate_ptr);
    const auto kPivotedLiteralPatternPtr 
              = pivoted_literal_pattern.FindVertex(gar_pattern_candidate_ptr->id());
    if (kPivotedLiteralPatternPtr) {
      // this vertex is contained in pivoted_literal_pattern
      auto [ pivoted_literal_pattern_candidate_set_it,
             pivoted_literal_pattern_candidate_set_ret ]
           = pivoted_literal_pattern_candidate_set.emplace(
            kPivotedLiteralPatternPtr,
            gar_pattern_candidate);
      // should added successfully
      assert(pivoted_literal_pattern_candidate_set_ret);
      const auto kPivotedRhsPatternPtr
                = pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id());
      if (kPivotedRhsPatternPtr){
        // this vertex is contained in pivoted_rhs_pattern
        auto [ pivoted_rhs_pattern_candidate_set_it,
               pivoted_rhs_pattern_candidate_set_ret ]
             = pivoted_rhs_pattern_candidate_set.emplace(
              kPivotedRhsPatternPtr,
                gar_pattern_candidate);
        // should added successfully
        assert(pivoted_rhs_pattern_candidate_set_ret);
        const auto kPivotedPatternPtr
                  = pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id());
        if (kPivotedPatternPtr){
          // this vertex is contained in pivoted_rhs_pattern
          auto [ pivoted_pattern_candidate_set_it,
                 pivoted_pattern_candidate_set_ret ]
               = pivoted_pattern_candidate_set.emplace(
                 kPivotedPatternPtr,
                  gar_pattern_candidate);
          // should added successfully
          assert(pivoted_pattern_candidate_set_ret);
        }
      }
      #ifndef NDEBUG
      else {
        // should not be contained in pivoted_rhs_pattern
        assert(!pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      }
      #endif // NDEBUG
    }
    #ifndef NDEBUG
    else {
      // should not be contained in pivoted_rhs_pattern and pivoted_pattern
      assert(!pivoted_rhs_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
      assert(    !pivoted_pattern.FindVertex(gar_pattern_candidate_ptr->id()));
    }
    #endif // NDEBUG
  }

  bool has_satisfy_x       = false,
       has_satisfy_x_not_y = false,
       has_satisfy_x_and_y = false;

  typename DataGraphType::VertexCounterType x_supp = 0,
                                           xy_supp = 0;



  /* ######################################################## *
   * ##  callbacks from pivoted_rhs_pattern to data graph  ## *
   * ######################################################## */
  // if the current match_state has already violate
  // a literal in x_literal_set, then prune it
  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_prune_callback
    = [&](const MatchType&
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    for (const auto& x_literal_ptr : pivoted_rhs_gar.x_literal_set()) {
      if (!x_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match)) {
        // this literal is not considered in the current
        // pivoted_rhs_pattern_to_data_graph_match, 
        // this match cannot be pruned move to the next literal
        continue;
      }
      if (!x_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        // does not satisfy all literals
        // does not need further matching
        return true;
      }
    }
    // satisfy all literals
    // continue matching
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_rhs_pattern_match_callback
    = [&](const MatchType& 
    pivoted_rhs_pattern_to_data_graph_match) -> bool {
    assert(!has_satisfy_x);
    assert(!has_satisfy_x_not_y);
    /* ############################################## *
     * ##  first check whether satisfy y literals  ## *
     * ############################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }

    // pivoted_rhs_gar_pattern is the same as gar_pattern
    has_satisfy_x = true;


    assert(!has_satisfy_x_not_y);
    assert(!has_satisfy_x_and_y);
    if (!has_satisfy_x) {
      // continue matching
      return true;
    }
    // satisfy x
    // reset parameter
    has_satisfy_x = false;
    if (satisify_y_literals) {
      // satisfy x and y
      assert(!has_satisfy_x_not_y);
      has_satisfy_x_and_y = true;
      // continue matching
      return true;
    }
    // satisfy x but not y
    has_satisfy_x_and_y = false;
    has_satisfy_x_not_y =  true;
    // no longer needs matching
    return false;
  };

  /* #################################################### *
   * ##  callbacks from pivoted_pattern to data graph  ## *
   * #################################################### */
  std::function<bool(const MatchType&)>
    pivoted_pattern_prune_callback
    = [](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    // prune nothing
    return false;
  };

  std::function<bool(const MatchType&)>
    pivoted_pattern_match_callback
    = [&](const MatchType& 
    pivoted_pattern_to_data_graph_match) -> bool {
    /* ######################################## *
     * ##  get a match for pivot vertex set  ## *
     * ######################################## */
    has_satisfy_x       = false;
    has_satisfy_x_not_y = false;
    has_satisfy_x_and_y = false;
    
    /* ####################################################### *
     * ##  get the paritial match from pivoted_rhs_pattern  ## *
     * ##  to data_graph                                    ## *
     * ####################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_rhs_pattern_to_data_graph_partial_match
      = pivoted_pattern_to_data_graph_match(
    pivoted_rhs_pattern_to_pivoted_pattern_partial_match);

    if (!pivoted_pattern_same_as_pivoted_rhs_pattern) {
      // pivoted_rhs_gar is not the same as literal_gar
      // first match the rhs literal, then match
      // the literal_gar
      std::cout << "should not be here" << std::endl;
      GUNDAM::MatchUsingMatch<match_semantics,
                              GUNDAM::MatchAlgorithm::kDagDp,
                              GUNDAM::MergeNecConfig::kNotMerge>(
                              pivoted_rhs_pattern,   data_graph, 
                              pivoted_rhs_pattern_to_data_graph_partial_match,
                              pivoted_rhs_pattern_candidate_set,
                              pivoted_rhs_pattern_prune_callback,
                              pivoted_rhs_pattern_match_callback,
                              time_limit_per_supp);
      assert(!has_satisfy_x);
      assert(!has_satisfy_x_and_y 
          || !has_satisfy_x_not_y); // cannot be both true
      has_satisfy_x = has_satisfy_x_and_y 
                   || has_satisfy_x_not_y;
      if (!has_satisfy_x) {
        return true;
      }
      x_supp++;
      const bool kSatisfyXSuppCallbackRet
                = satisfy_x_supp_callback(pivoted_pattern_to_data_graph_match);
      if (has_satisfy_x_not_y) {
        return kSatisfyXSuppCallbackRet;
      }
      assert(has_satisfy_x_and_y);
      xy_supp++;
      const bool kSatisfyXYSuppCallbackRet
                = satisfy_xy_supp_callback(pivoted_pattern_to_data_graph_match);
      return kSatisfyXSuppCallbackRet
         && kSatisfyXYSuppCallbackRet;
    }
    /* ########################################################### *
     * ##  pivoted_pattern is the same as pivoted_rhs_pattern,  ## *
     * ##  directly whether the match satisfy y literals here   ## *
     * ########################################################### */
    assert(pivoted_rhs_pattern_to_data_graph_partial_match.size()
            == pivoted_pattern_to_data_graph_match.size());

    /* ########################################################### *
     * ##  get the paritial match from pivoted_literal_pattern  ## *
     * ##  to data_graph                                        ## *
     * ########################################################### */
    GUNDAM::Match<GraphPatternType, DataGraphType>
    pivoted_literal_pattern_to_data_graph_partial_match
          = pivoted_pattern_to_data_graph_match(
    pivoted_literal_pattern_to_pivoted_pattern_partial_match);

    // pivoted_rhs_pattern is the same as pivoted_literal_pattern
    // directly match to pivoted_literal_pattern
    /* ############################################### *
      * ##  get the paritial match from gar_pattern  ## *
      * ##  to data_graph                            ## *
      * ############################################### */
    assert(        SameGar    (pivoted_rhs_gar, pivoted_literal_gar));
    assert(GUNDAM::SamePattern(pivoted_pattern, pivoted_literal_pattern));

    /* ######################################################## *
    * ##  chech whether the match satisfy rhs literal here  ## *
    * ######################################################## */
    bool satisify_y_literals = true;
    for (const auto& y_literal_ptr : pivoted_rhs_gar.y_literal_set()) {
      assert(y_literal_ptr->MappedBy(pivoted_rhs_pattern_to_data_graph_partial_match));
      if (y_literal_ptr->Satisfy(pivoted_rhs_pattern_to_data_graph_partial_match)) {
        continue;
      }
      satisify_y_literals = false;
      break;
    }
    if (!satisify_y_literals) {
      return true;
    }
    // satisfy y literals
    xy_supp++;
    // continue matching
    const bool kSatisfyXYSuppCallbackRet
              = satisfy_xy_supp_callback(gar_pattern_to_data_graph_partial_match);
    return true;
  };

  MatchType pivoted_pattern_partial_match;
  for (auto pivoted_pattern_vertex_it = pivoted_pattern.VertexBegin();
           !pivoted_pattern_vertex_it.IsDone();
            pivoted_pattern_vertex_it++) {
    auto target_handle = gar_pattern_to_data_graph_partial_match.MapTo(pivoted_pattern_vertex_it->id());
    if (!target_handle) {
      continue;
    }
    pivoted_pattern_partial_match.AddMap(pivoted_pattern_vertex_it,
                                         target_handle);
  }

  // first match pivoted gar_pattern, then match rhs gar_pattern
  GUNDAM::MatchUsingMatch<match_semantics,
                          GUNDAM::MatchAlgorithm::kDagDp,
                          GUNDAM::MergeNecConfig::kNotMerge>(
                          pivoted_pattern, data_graph,
                          pivoted_pattern_partial_match,
                          pivoted_pattern_candidate_set,
                          pivoted_pattern_prune_callback,
                          pivoted_pattern_match_callback,
                          time_limit);
  assert(x_supp >= xy_supp);
  return std::pair(x_supp, xy_supp);
}


double ComputeProb(double p, int64_t t, int64_t trial_num) {
  double log_prob = 0.00;
  if (t == 0) {
    log_prob += log(1 - p) * trial_num;
    return exp(log_prob);
  }
  if (t == trial_num) {
    log_prob += log(p) * trial_num;
    return exp(log_prob);
  }
  for (int64_t i = trial_num - t + 1; i <= trial_num; i++) {
    log_prob += log(i);
  }
  for (int64_t i = 1; i <= t; i++) {
    log_prob -= log(i);
  }
  log_prob += log(p) * t;
  log_prob += log(1 - p) * (trial_num - t);
  return exp(log_prob);
}



int64_t EstimateUpperBound(double success_probability, int64_t trial_num, double threshold) {
  double fail_probability = 1.0 - success_probability;
  double current_prob = 0.0;
  for (int64_t t = 0; t <= trial_num; t++) {
    double prob = ComputeProb(fail_probability, t, trial_num);
    if (current_prob + prob >= (1.0 - threshold )) {
      return trial_num - t;
    }
  }
  return 0;
}

int64_t EstimateLowerBound(double success_probability, int64_t trial_num, double threshold) {
  double current_prob = 0.0;
  for (int64_t t = 0; t <= trial_num; t++) {
    double prob = ComputeProb(success_probability, t, trial_num);
    if (current_prob + prob >= (1.0 - threshold )) {
      return t;
    }
  }
  return trial_num;
}



} // namespace gar

#endif // GAR_SUPP_H_
