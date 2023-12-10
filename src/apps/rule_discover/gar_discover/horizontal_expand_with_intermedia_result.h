#ifndef EXAMPLES_ANALYTICAL_HORIZONTAL_EXPAND_WITH_INTERMEDIA_RESULT_
#define EXAMPLES_ANALYTICAL_HORIZONTAL_EXPAND_WITH_INTERMEDIA_RESULT_

#include "generate_tree.h"

#include "gar/literal_info.h"
#include "gar/literal_stand_alone_info.h"

template <typename GraphPatternType,
          typename    DataGraphType, 
inline void HorizontalExpandWithIntermediaResult (
     const std::vector<
      gar::LiteralInfo<GraphPatternType,
                          DataGraphType>>& rhs_literal_info_set,
     const std::vector<
      gar::LiteralInfo<GraphPatternType,
                          DataGraphType>>& lhs_literal_info_not_in_rhs_set,
     const std::vector<
      std::map<typename GUNDAM::VertexHandle<GraphPatternType>::type,
      std::vector<typename GUNDAM::VertexHandle<DataGraphType>::type>>>& candidate_set_for_data_graph_set,
        ExpandTreeNode<GraphPatternType,
                          DataGraphType>& expand_node,
      GenerateTreeNode<GraphPatternType,
                          DataGraphType>& generate_node ) {
  for (const auto& rhs_literal 
                 : rhs_literal_info_set) {
    // do something here
    generate_node.AddLiteralTree(rhs_literal,
                                 rhs_literal_conf);
  }
  return;
}

#endif // EXAMPLES_ANALYTICAL_HORIZONTAL_EXPAND_WITH_INTERMEDIA_RESULT_