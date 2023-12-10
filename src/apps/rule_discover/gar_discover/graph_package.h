#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GRAPH_PACKAGE_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GRAPH_PACKAGE_H_

#include "gundam/algorithm/neighborhood_equivalence_class.h"

#include "gundam/tool/isolate_vertex.h"
#include "gundam/tool/vertex_degree_filter.h"
#include "gundam/data_type/datatype.h"

#include "gundam/io/csvgraph.h"

#include "gundam/graph_statistics/graph_basic_statistics.h"

#include "gundam/type_getter/vertex_label.h"
#include "gundam/type_getter/vertex_id.h"
#include "gundam/type_getter/edge_label.h"
#include "gundam/type_getter/edge_id.h"

#include "util/log.h"

namespace grape {

namespace _gar_discover { 

template <typename GraphType>
class GraphPackage{
 private:
  using VertexIDType = typename GUNDAM::VertexID<GraphType>::type;
  using   EdgeIDType = typename GUNDAM::  EdgeID<GraphType>::type;

  using VertexLabelType = typename GUNDAM::VertexLabel<GraphType>::type;
  using   EdgeLabelType = typename GUNDAM::  EdgeLabel<GraphType>::type;

  using EdgeTypeType = std::tuple<VertexLabelType,
                                    EdgeLabelType,
                                  VertexLabelType>;

  using GraphBasicStatisticsType = GUNDAM::GraphBasicStatistics<GraphType>;

 public:
  GraphPackage(
    const std::set<EdgeTypeType>&  specified_edge_type_set,
    const std::set<EdgeLabelType>& specified_edge_label_set,
    const double constant_freq_bound,
    const std::string& kGraphPathVFile,
    const std::string& kGraphPathEFile,
    GraphBasicStatisticsType& graph_basic_statistics,
    const std::string& kKnowledgeGraphPathVFile = "",
    const std::string& kKnowledgeGraphPathEFile = "",
    const std::string& kERFile = "") {

    this->load_data_graph_success_ = false;
      
    this->graph_path_v_file_ = kGraphPathVFile;
    this->graph_path_e_file_ = kGraphPathEFile;


    if (kKnowledgeGraphPathVFile == ""
        && kKnowledgeGraphPathEFile == "") {
      if (GUNDAM::ReadCSVGraph(this->data_graph_, 
                                kGraphPathVFile,
                                kGraphPathEFile) < 0) {
        util::Error("load data graph failed!");
        return;
      }
    } else {
      if (GUNDAM::ReadTwoCSVGraph(this->data_graph_, 
                                kGraphPathVFile,
                                kGraphPathEFile,
                                kKnowledgeGraphPathVFile,
                                kKnowledgeGraphPathEFile,
                                kERFile) < 0) {
        util::Error("load data and knowledge graph failed");
        return;
      }
    }

    this->load_data_graph_success_ = true;

    graph_basic_statistics.AddGraph(this->data_graph_);

    const auto kMaxVertexId = graph_basic_statistics.max_vertex_id();
    const auto kMaxEdgeId   = graph_basic_statistics.max_edge_id();

    util::Info(" vertex_label_set.size(): " + std::to_string(graph_basic_statistics.vertex_label_counter().size()));
    util::Info("   edge_label_set.size(): " + std::to_string(graph_basic_statistics.  edge_label_counter().size()));
    util::Info("    edge_type_set.size(): " + std::to_string(graph_basic_statistics.   edge_type_counter().size()));

    this->max_edge_id_ = kMaxEdgeId;

    if (!specified_edge_type_set.empty()) {

      graph_basic_statistics.PreserveEdgeTypeSet(specified_edge_type_set);

      util::Info("################################");
      util::Info("##   preserve edge type set   ##");
      util::Info("################################");

      util::Info("   edge type considered: "
                + std::to_string(graph_basic_statistics.edge_type_counter().size()));
      util::Info("vertex label considered: "
                + std::to_string(graph_basic_statistics.vertex_label_counter().size()));
      util::Info("  edge label considered: "
                + std::to_string(graph_basic_statistics.edge_label_counter().size()));
    }

    if (!specified_edge_label_set.empty()) {
      graph_basic_statistics.PreserveEdgeLabelSet(specified_edge_label_set);

      util::Info("#################################");
      util::Info("##   preserve edge label set   ##");
      util::Info("#################################");

      util::Info("   edge type considered: "
                + std::to_string(graph_basic_statistics.edge_type_counter().size()));
      util::Info("vertex label considered: "
                + std::to_string(graph_basic_statistics.vertex_label_counter().size()));
      util::Info("  edge label considered: "
                + std::to_string(graph_basic_statistics.edge_label_counter().size()));
    }
    
    if (!specified_edge_type_set.empty()
     || !specified_edge_label_set.empty()) {
        
      if constexpr (GUNDAM::GraphParameter<GraphType>::graph_level_erase_vertex_iterator) {
        for (auto vertex_it = this->data_graph_.VertexBegin();
                 !vertex_it.IsDone();) {
          if (graph_basic_statistics.vertex_label_counter().find(vertex_it->label())
           == graph_basic_statistics.vertex_label_counter().end()) {
            // erase
            vertex_it = this->data_graph_.EraseVertex(vertex_it);
            continue;
          }
          // preserve
          vertex_it++;
        }
      }
      else {
        // cannot call:
        //    vertex_it = this->data_graph_.EraseVertex(vertex_it);
        std::vector<VertexIDType> vertex_id_to_remove;
        vertex_id_to_remove.reserve(this->data_graph_.CountVertex());
        // just begin bfs at a random vertex and find whether
        // it can reach all vertexes
        for (auto vertex_it = this->data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++){
          if (graph_basic_statistics.vertex_label_counter().find(vertex_it->label())
           != graph_basic_statistics.vertex_label_counter().end()) {
            // preserve
            continue;
          }
          // erase
          vertex_id_to_remove.emplace_back(vertex_it->id());
        }
        #ifndef NDEBUG 
        const size_t kVertexNum = this->data_graph_.CountVertex();
        #endif // NDEBUG
        for (const auto& vertex_id : vertex_id_to_remove) {
          this->data_graph_.EraseVertex(vertex_id);
        }
        assert(kVertexNum == this->data_graph_.CountVertex() + vertex_id_to_remove.size());
      }

      if constexpr (GUNDAM::GraphParameter<GraphType>::vertex_level_erase_edge_iterator) {
        for (auto vertex_it = this->data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          for (auto out_edge_it = vertex_it->OutEdgeBegin();
                   !out_edge_it.IsDone();) {
            const auto edge_type = std::tuple(out_edge_it->src_handle()->label(),
                                              out_edge_it->label(),
                                              out_edge_it->dst_handle()->label());
            if (graph_basic_statistics.edge_type_counter().find(edge_type)
             != graph_basic_statistics.edge_type_counter().end()) {
              // this edge does not need to be removed;
              out_edge_it++;
              continue;
            }
            out_edge_it = vertex_it->EraseEdge(out_edge_it);
          }
        }
      }
      else {
        // cannot call:
        //    out_edge_it = vertex_it->EraseEdge(out_edge_it);
        std::vector<EdgeIDType> edge_id_to_remove;
        edge_id_to_remove.reserve(this->data_graph_.CountEdge());
        for (auto vertex_it = this->data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          for (auto out_edge_it = vertex_it->OutEdgeBegin();
                   !out_edge_it.IsDone();) {
            const auto edge_type = std::tuple(out_edge_it->src_handle()->label(),
                                              out_edge_it->label(),
                                              out_edge_it->dst_handle()->label());
            if (graph_basic_statistics.edge_type_counter().find(edge_type)
             != graph_basic_statistics.edge_type_counter().end()) {
              // this edge does not need to be removed;
              out_edge_it++;
              continue;
            }
            edge_id_to_remove.emplace_back(out_edge_it->id());
          }
        }
        #ifndef NDEBUG 
        const size_t kEdgeNum = this->data_graph_.CountEdge();
        #endif // NDEBUG
        for (const auto& edge_id : edge_id_to_remove) {
          this->data_graph_.EraseEdge(edge_id);
        }
        assert(kEdgeNum == this->data_graph_.CountEdge() + edge_id_to_remove.size());
      }

      GUNDAM::RemoveIsolateVertex(this->data_graph_);
    }

    // assume that the vertex with the same label has the same
    // attributes set, collect all legal keys and values for
    // attribute
    graph_basic_statistics.ResetVertexAttr(this->data_graph_,
                                            constant_freq_bound);

    util::Info("data graph vertex number: " + std::to_string(this->data_graph_.CountVertex()));
    util::Info("data graph  edge  number: " + std::to_string(this->data_graph_.CountEdge  ()));
    return;
  }

  inline const std::string& graph_path_v_file() const {
    return this->graph_path_v_file_;
  }

  inline const std::string& graph_path_e_file() const {
    return this->graph_path_e_file_;
  }

  inline bool AddMlEdge(const std::string& kMlLiteralEdgesFile,
                   std::map<EdgeLabelType,
                            EdgeLabelType>& normal_to_ml_edge_label) {
    if (kMlLiteralEdgesFile == "") {
      // does not have ml edge
      return true;
    }
    // has ml literal
    std::ifstream ml_literal_edge_file(kMlLiteralEdgesFile);
    
    if (!ml_literal_edge_file.good()){
      std::cout << " ml literal edge file is not good! " << std::endl;
      return false;
    }
                  
    EdgeIDType edge_id_allocator = this->max_edge_id_;
    edge_id_allocator++;

    while (ml_literal_edge_file) {
      std::string s;
      if (!std::getline( ml_literal_edge_file, s )) 
        break;

      std::istringstream ss( s );
      std::vector <std::string> record;
      record.reserve(3);
      while (ss) {
        std::string s;
        if (!std::getline( ss, s, ',' )) {
          break;
        }
        record.emplace_back( s );
      }
       VertexIDType  src_id    = GUNDAM::StringToDataType< VertexIDType>(record[0]);
      EdgeLabelType edge_label = GUNDAM::StringToDataType<EdgeLabelType>(record[1]);
       VertexIDType  dst_id    = GUNDAM::StringToDataType< VertexIDType>(record[2]);

      assert(normal_to_ml_edge_label.find(edge_label)
          != normal_to_ml_edge_label.end());
      EdgeLabelType ml_edge_label = normal_to_ml_edge_label[edge_label];
      auto [ edge_handle,
              edge_ret ] = this->data_graph_.AddEdge(src_id, dst_id, 
                                                  ml_edge_label, 
                                                    edge_id_allocator++);
      assert(edge_ret);
    }
    return true;
  }

  GraphType& data_graph(){
    return this->data_graph_;
  }

  inline const auto& data_graph_nec() const {
    return this->data_graph_nec_;
  }

  inline void GenerateGraphNec() {
    this->data_graph_nec_ = GUNDAM::Nec(this->data_graph_);
    return;
  }

  inline bool load_data_graph_success() const {
    return this->load_data_graph_success_;
  }

  inline void GenerateCentralSet() {
    std::vector<typename GUNDAM::VertexHandle<GraphType>::type> centrel_vertex;
    centrel_vertex = GUNDAM::VertexDegreeFilter<
                     GUNDAM::FilterType::kHigherOrEqualTo,
                     GUNDAM::EdgeDirection::kInOut>(this->data_graph_, 3);
    for (const auto& vertex : centrel_vertex) {
      this->centrel_vertex_set_[vertex->label()].emplace_back(vertex);
    }
    for (auto& [vertex_label, centrel_vertex_set] 
                      : this->centrel_vertex_set_) {
      std::sort(centrel_vertex_set.begin(),
                centrel_vertex_set.end());
      util::Debug("vertex_label: "              + std::to_string(vertex_label));
      util::Debug("centrel_vertex_set.size(): " + std::to_string(centrel_vertex_set.size()));
    }
    return;
  }

  inline const auto& centrel_vertex_set(
         const typename GUNDAM::VertexLabel <GraphType>::type& vertex_label) const {
    auto centrel_vertex_set_it = this->centrel_vertex_set_.find(vertex_label);
    assert(centrel_vertex_set_it
        != this->centrel_vertex_set_.end());
    return centrel_vertex_set_it->second;
  }

  inline void SampleCentralSet(double sample_ratio) {
    assert(!this->centrel_vertex_set_.empty());
    assert(sample_ratio > 0 
        && sample_ratio < 1.0);
    for (auto& [vertex_label, centrel_vertex_set] 
                      : this->centrel_vertex_set_) {
      assert(std::is_sorted(centrel_vertex_set.begin(),
                            centrel_vertex_set.end()));
      std::vector<
      typename GUNDAM::VertexHandle<GraphType>::type> preserved_centrel_vertex_set;
      for (const auto& centrel_vertex_handle  
                     : centrel_vertex_set) {
        const double kR = ((double) rand() / (RAND_MAX));
        if (kR > sample_ratio) {
          /// this vertex should not be preserved
          continue;
        }
        /// this vertex should be preserved
        preserved_centrel_vertex_set.emplace_back(centrel_vertex_handle);
      }
      centrel_vertex_set.swap(preserved_centrel_vertex_set);
      assert(std::is_sorted(centrel_vertex_set.begin(),
                            centrel_vertex_set.end()));
      util::Debug("vertex_label: "              + std::to_string(vertex_label));
      util::Debug("centrel_vertex_set.size(): " + std::to_string(centrel_vertex_set.size()));
    }
    return;
  }

 private:
  bool load_data_graph_success_;

  std::string graph_path_v_file_,
              graph_path_e_file_;

   GraphType  data_graph_;
  EdgeIDType max_edge_id_;

  std::vector<
  std::vector<typename GUNDAM::VertexHandle<GraphType>::type>> data_graph_nec_;

  std::map<typename GUNDAM::VertexLabel <GraphType>::type,
           std::vector<
           typename GUNDAM::VertexHandle<GraphType>::type>> centrel_vertex_set_;
};



template <typename GraphType>
class IncGraphPackage{
 private:
  using VertexIDType = typename GUNDAM::VertexID<GraphType>::type;
  using   EdgeIDType = typename GUNDAM::  EdgeID<GraphType>::type;

  using VertexLabelType = typename GUNDAM::VertexLabel<GraphType>::type;
  using   EdgeLabelType = typename GUNDAM::  EdgeLabel<GraphType>::type;

  using EdgeTypeType = std::tuple<VertexLabelType,
                                    EdgeLabelType,
                                  VertexLabelType>;

  using GraphBasicStatisticsType = GUNDAM::GraphBasicStatistics<GraphType>;

 public:
  IncGraphPackage(
    const std::set<EdgeTypeType>&  specified_edge_type_set,
    const std::set<EdgeLabelType>& specified_edge_label_set,
    const double constant_freq_bound,
    const std::string& kGraphPathVFile,
    const std::string& kGraphPathEFile,
    // data graph G
    GraphBasicStatisticsType& origin_graph_basic_statistics, 
    // data graph G \oplus \Delta G
    GraphBasicStatisticsType& updated_graph_basic_statistics,
    const std::string& kKnowledgeGraphVFile = "",
    const std::string& kKnowledgeGraphEFile = "",
    const std::string& ERFile = "") {

    this->load_data_graph_success_ = false;
      
    this->graph_path_v_file_ = kGraphPathVFile;
    this->graph_path_e_file_ = kGraphPathEFile;

    if ((kKnowledgeGraphVFile == "")
        && (kKnowledgeGraphEFile == "")) {
      if (GUNDAM::ReadIncCSVGraph(this->origin_data_graph_,
                                this->updated_data_graph_,
                                this->id_to_origin_updated_vertex_handle_,
                                this->updated_vertex_id_,
                                kGraphPathVFile,
                                kGraphPathEFile) < 0) {
        util::Error("load incremental data graph failed!");
        return;
      }
    } else {
      if (GUNDAM::ReadTwoIncCSVGraph(this->origin_data_graph_,
                                this->updated_data_graph_,
                                this->id_to_origin_updated_vertex_handle_,
                                this->updated_vertex_id_,
                                kGraphPathVFile,
                                kGraphPathEFile,
                                kKnowledgeGraphVFile,
                                kKnowledgeGraphEFile,
                                ERFile) < 0) {
        util::Error("load two incremental data graph failed!");
        return;
      }
    }

    this->load_data_graph_success_ = true;

    origin_graph_basic_statistics.AddGraph(this->origin_data_graph_);
    updated_graph_basic_statistics.AddGraph(this->updated_data_graph_);

    const auto kOriginMaxVertexId = origin_graph_basic_statistics.max_vertex_id();
    const auto kOriginMaxEdgeId   = origin_graph_basic_statistics.max_edge_id();

    const auto kUpdatedMaxVertexId = updated_graph_basic_statistics.max_vertex_id();
    const auto kUpdatedMaxEdgeId   = updated_graph_basic_statistics.max_edge_id();

    // those are same in origin graph and updated graph

    util::Info(" vertex_label_set.size(): " + std::to_string(origin_graph_basic_statistics.vertex_label_counter().size()));
    util::Info("   edge_label_set.size(): " + std::to_string(origin_graph_basic_statistics.  edge_label_counter().size()));
    util::Info("    edge_type_set.size(): " + std::to_string(origin_graph_basic_statistics.   edge_type_counter().size()));

    this->origin_max_edge_id_ = kOriginMaxEdgeId;
    this->updated_max_edge_id_ = kUpdatedMaxEdgeId;

    if (!specified_edge_type_set.empty()) {

      origin_graph_basic_statistics.PreserveEdgeTypeSet(specified_edge_type_set);
      updated_graph_basic_statistics.PreserveEdgeTypeSet(specified_edge_type_set);

      util::Info("################################");
      util::Info("##   preserve edge type set   ##");
      util::Info("################################");

      util::Info("   edge type considered in origin graph: "
                + std::to_string(origin_graph_basic_statistics.edge_type_counter().size()));
      util::Info("vertex label considered in origin graph: "
                + std::to_string(origin_graph_basic_statistics.vertex_label_counter().size()));
      util::Info("  edge label considered in origin graph: "
                + std::to_string(origin_graph_basic_statistics.edge_label_counter().size()));

      util::Info("   edge type considered in updated graph: "
                + std::to_string(updated_graph_basic_statistics.edge_type_counter().size()));
      util::Info("vertex label considered in updated graph: "
                + std::to_string(updated_graph_basic_statistics.vertex_label_counter().size()));
      util::Info("  edge label considered in updated graph: "
                + std::to_string(updated_graph_basic_statistics.edge_label_counter().size()));
    }

    if (!specified_edge_label_set.empty()) {
      origin_graph_basic_statistics.PreserveEdgeLabelSet(specified_edge_label_set);
      updated_graph_basic_statistics.PreserveEdgeLabelSet(specified_edge_label_set);

      util::Info("#################################");
      util::Info("##   preserve edge label set   ##");
      util::Info("#################################");

      util::Info("  origin edge type considered: "
                + std::to_string(origin_graph_basic_statistics.edge_type_counter().size()));
      util::Info("  origin vertex label considered: "
                + std::to_string(origin_graph_basic_statistics.vertex_label_counter().size()));
      util::Info("  origin edge label considered: "
                + std::to_string(origin_graph_basic_statistics.edge_label_counter().size()));

      util::Info("  updated edge type considered: "
                + std::to_string(updated_graph_basic_statistics.edge_type_counter().size()));
      util::Info("  updated vertex label considered: "
                + std::to_string(updated_graph_basic_statistics.vertex_label_counter().size()));
      util::Info("  updated edge label considered: "
                + std::to_string(updated_graph_basic_statistics.edge_label_counter().size()));
    }
    
    if (!specified_edge_type_set.empty()
     || !specified_edge_label_set.empty()) {
      std::vector<VertexIDType> origin_remove_vertex_id, updated_remove_vertex_id;
      if constexpr (GUNDAM::GraphParameter<GraphType>::graph_level_erase_vertex_iterator) {
        for (auto vertex_it = this->origin_data_graph_.VertexBegin();
                 !vertex_it.IsDone();) {
          if (origin_graph_basic_statistics.vertex_label_counter().find(vertex_it->label())
           == origin_graph_basic_statistics.vertex_label_counter().end()) {
            // erase
            vertex_it = this->origin_data_graph_.EraseVertex(vertex_it);
            origin_remove_vertex_id.emplace_back(vertex_it->id());
            continue;
          }
          // preserve
          vertex_it++;
        }

        for (auto vertex_it = this->updated_data_graph_.VertexBegin();
                 !vertex_it.IsDone();) {
          if (updated_graph_basic_statistics.vertex_label_counter().find(vertex_it->label())
           == updated_graph_basic_statistics.vertex_label_counter().end()) {
            // erase
            vertex_it = this->updated_data_graph_.EraseVertex(vertex_it);
            updated_remove_vertex_id.emplace_back(vertex_it->id());
            continue;
          }
          // preserve
          vertex_it++;
        }
      } else {
        // cannot call:
        //    vertex_it = this->data_graph_.EraseVertex(vertex_it);
        std::vector<VertexIDType> vertex_id_to_remove;
        vertex_id_to_remove.reserve(this->origin_data_graph_.CountVertex());
        // just begin bfs at a random vertex and find whether
        // it can reach all vertexes
        for (auto vertex_it = this->origin_data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          if (origin_graph_basic_statistics.vertex_label_counter().find(vertex_it->label())
           != origin_graph_basic_statistics.vertex_label_counter().end()) {
            // preserve
            continue;
          }
          // erase
          vertex_id_to_remove.emplace_back(vertex_it->id());
          origin_remove_vertex_id.emplace_back(vertex_it->id());
        }
        #ifndef NDEBUG 
        const size_t kOriginVertexNum = this->origin_data_graph_.CountVertex();
        #endif // NDEBUG
        for (const auto& vertex_id : vertex_id_to_remove) {
          this->origin_data_graph_.EraseVertex(vertex_id);
        }
        assert(kOriginVertexNum == this->origin_data_graph_.CountVertex() + vertex_id_to_remove.size());

        vertex_id_to_remove.clear();
        vertex_id_to_remove.reserve(this->updated_data_graph_.CountVertex());
        // just begin bfs at a random vertex and find whether
        // it can reach all vertexes
        for (auto vertex_it = this->updated_data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          if (updated_graph_basic_statistics.vertex_label_counter().find(vertex_it->label())
           != updated_graph_basic_statistics.vertex_label_counter().end()) {
            // preserve
            continue;
          }
          // erase
          vertex_id_to_remove.emplace_back(vertex_it->id());
          updated_remove_vertex_id.emplace_back(vertex_it->id());
        }
        #ifndef NDEBUG 
        const size_t kUpdatedVertexNum = this->updated_data_graph_.CountVertex();
        #endif // NDEBUG
        for (const auto& vertex_id : vertex_id_to_remove) {
          this->updated_data_graph_.EraseVertex(vertex_id);
        }
        assert(kUpdatedVertexNum == this->updated_data_graph_.CountVertex() + vertex_id_to_remove.size());
      }

      for (auto vertex_id : origin_remove_vertex_id) {
        auto updated_handle = id_to_origin_updated_vertex_handle_[vertex_id].second;
        id_to_origin_updated_vertex_handle_[vertex_id] = std::make_pair(nullptr, updated_handle);
      }

      for (auto vertex_id : updated_remove_vertex_id) {
        auto origin_handle = id_to_origin_updated_vertex_handle_[vertex_id].first;
        id_to_origin_updated_vertex_handle_[vertex_id] = std::make_pair(origin_handle, nullptr);
      }

      if constexpr (GUNDAM::GraphParameter<GraphType>::vertex_level_erase_edge_iterator) {
        for (auto vertex_it = this->origin_data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          for (auto out_edge_it = vertex_it->OutEdgeBegin();
                   !out_edge_it.IsDone();) {
            const auto edge_type = std::tuple(out_edge_it->src_handle()->label(),
                                              out_edge_it->label(),
                                              out_edge_it->dst_handle()->label());
            if (origin_graph_basic_statistics.edge_type_counter().find(edge_type)
             != origin_graph_basic_statistics.edge_type_counter().end()) {
              // this edge does not need to be removed;
              out_edge_it++;
              continue;
            }
            out_edge_it = vertex_it->EraseEdge(out_edge_it);
          }
        }

        for (auto vertex_it = this->updated_data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          for (auto out_edge_it = vertex_it->OutEdgeBegin();
                   !out_edge_it.IsDone();) {
            const auto edge_type = std::tuple(out_edge_it->src_handle()->label(),
                                              out_edge_it->label(),
                                              out_edge_it->dst_handle()->label());
            if (updated_graph_basic_statistics.edge_type_counter().find(edge_type)
             != updated_graph_basic_statistics.edge_type_counter().end()) {
              // this edge does not need to be removed;
              out_edge_it++;
              continue;
            }
            out_edge_it = vertex_it->EraseEdge(out_edge_it);
          }
        }        
      } else {
        // cannot call:
        //    out_edge_it = vertex_it->EraseEdge(out_edge_it);
        std::vector<EdgeIDType> edge_id_to_remove;
        edge_id_to_remove.reserve(this->origin_data_graph_.CountEdge());
        for (auto vertex_it = this->origin_data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          for (auto out_edge_it = vertex_it->OutEdgeBegin();
                   !out_edge_it.IsDone();) {
            const auto edge_type = std::tuple(out_edge_it->src_handle()->label(),
                                              out_edge_it->label(),
                                              out_edge_it->dst_handle()->label());
            if (origin_graph_basic_statistics.edge_type_counter().find(edge_type)
             != origin_graph_basic_statistics.edge_type_counter().end()) {
              // this edge does not need to be removed;
              out_edge_it++;
              continue;
            }
            edge_id_to_remove.emplace_back(out_edge_it->id());
          }
        }
        #ifndef NDEBUG 
        const size_t kOriginEdgeNum = this->origin_data_graph_.CountEdge();
        #endif // NDEBUG
        for (const auto& edge_id : edge_id_to_remove) {
          this->origin_data_graph_.EraseEdge(edge_id);
        }
        assert(kOriginEdgeNum == this->origin_data_graph_.CountEdge() + edge_id_to_remove.size());

        edge_id_to_remove.reserve(this->updated_data_graph_.CountEdge());
        for (auto vertex_it = this->updated_data_graph_.VertexBegin();
                 !vertex_it.IsDone();
                  vertex_it++) {
          for (auto out_edge_it = vertex_it->OutEdgeBegin();
                   !out_edge_it.IsDone();) {
            const auto edge_type = std::tuple(out_edge_it->src_handle()->label(),
                                              out_edge_it->label(),
                                              out_edge_it->dst_handle()->label());
            if (updated_graph_basic_statistics.edge_type_counter().find(edge_type)
             != updated_graph_basic_statistics.edge_type_counter().end()) {
              // this edge does not need to be removed;
              out_edge_it++;
              continue;
            }
            edge_id_to_remove.emplace_back(out_edge_it->id());
          }
        }
        #ifndef NDEBUG 
        const size_t kUpdatedEdgeNum = this->updated_data_graph_.CountEdge();
        #endif // NDEBUG
        for (const auto& edge_id : edge_id_to_remove) {
          this->updated_data_graph_.EraseEdge(edge_id);
        }
        assert(kUpdatedEdgeNum == this->updated_data_graph_.CountEdge() + edge_id_to_remove.size());        
      }

      std::set<VertexIDType> removed_id_set;
      GUNDAM::RemoveIsolateVertex(this->origin_data_graph_, removed_id_set);

      for (auto vertex_id : removed_id_set) {
        auto updated_handle = id_to_origin_updated_vertex_handle_[vertex_id].second;
        id_to_origin_updated_vertex_handle_[vertex_id] = std::make_pair(nullptr, updated_handle);
      }

      removed_id_set.clear();
      GUNDAM::RemoveIsolateVertex(this->updated_data_graph_, removed_id_set);

      for (auto vertex_id : removed_id_set) {
        auto origin_handle = id_to_origin_updated_vertex_handle_[vertex_id].first;
        id_to_origin_updated_vertex_handle_[vertex_id] = std::make_pair(origin_handle, nullptr);
      }
    }

    // assume that the vertex with the same label has the same
    // attributes set, collect all legal keys and values for
    // attribute
    origin_graph_basic_statistics.ResetVertexAttr(this->origin_data_graph_,
                                            constant_freq_bound);
    updated_graph_basic_statistics.ResetVertexAttr(this->updated_data_graph_,
                                            constant_freq_bound);

    util::Info("origin graph vertex number: " + std::to_string(this->origin_data_graph_.CountVertex()));
    util::Info("origin graph  edge  number: " + std::to_string(this->origin_data_graph_.CountEdge  ()));

    util::Info("updated graph vertex number: " + std::to_string(this->updated_data_graph_.CountVertex()));
    util::Info("updated graph  edge  number: " + std::to_string(this->updated_data_graph_.CountEdge  ()));    
    return;
  }

  inline const std::string& graph_path_v_file() const {
    return this->graph_path_v_file_;
  }

  inline const std::string& graph_path_e_file() const {
    return this->graph_path_e_file_;
  }

  inline bool AddMlEdge(const std::string& kMlLiteralEdgesFile,
                   std::map<EdgeLabelType,
                            EdgeLabelType>& normal_to_ml_edge_label) {
    if (kMlLiteralEdgesFile == "") {
      // does not have ml edge
      return true;
    }
    // has ml literal
    std::ifstream ml_literal_edge_file(kMlLiteralEdgesFile);
    
    if (!ml_literal_edge_file.good()){
      std::cout << " ml literal edge file is not good! " << std::endl;
      return false;
    }
                  
    EdgeIDType origin_edge_id_allocator = this->origin_max_edge_id_,
              updated_edge_id_allocator = this->updated_max_edge_id_;
    origin_edge_id_allocator++;
    updated_edge_id_allocator++;

    while (ml_literal_edge_file) {
      std::string s;
      if (!std::getline( ml_literal_edge_file, s )) 
        break;

      std::istringstream ss( s );
      std::vector <std::string> record;
      record.reserve(3);
      while (ss) {
        std::string s;
        if (!std::getline( ss, s, ',' )) {
          break;
        }
        record.emplace_back( s );
      }
       VertexIDType  src_id    = GUNDAM::StringToDataType< VertexIDType>(record[0]);
      EdgeLabelType edge_label = GUNDAM::StringToDataType<EdgeLabelType>(record[1]);
       VertexIDType  dst_id    = GUNDAM::StringToDataType< VertexIDType>(record[2]);


      // assume the ml edge is not affected by inc
      assert(normal_to_ml_edge_label.find(edge_label)
          != normal_to_ml_edge_label.end());
      EdgeLabelType ml_edge_label = normal_to_ml_edge_label[edge_label];
      auto [ edge_handle,
              edge_ret ] = this->origin_data_graph_.AddEdge(src_id, dst_id, 
                                                  ml_edge_label, 
                                                    origin_edge_id_allocator++);
      assert(edge_ret);
      auto [updated_edge_handle,
              updated_edge_ret] = this->updated_data_graph_.AddEdge(src_id, dst_id, 
                                                  ml_edge_label, 
                                                    updated_edge_id_allocator++);
      assert(updated_edge_ret);
    }
    return true;
  }

  GraphType& origin_data_graph(){
    return this->origin_data_graph_;
  }

  GraphType& updated_data_graph(){
    return this->updated_data_graph_;
  }

  inline const auto& origin_data_graph_nec() const {
    return this->origin_data_graph_nec_;
  }

  inline const auto& updated_data_graph_nec() const {
    return this->updated_data_graph_nec_;
  }

  inline const auto& id_to_vertex_handle() const {
    return this->id_to_origin_updated_vertex_handle_;
  }

  inline const auto& updated_vertex_id() const {
    return this->updated_vertex_id_;
  }

  inline void GenerateGraphNec() {
    this->origin_data_graph_nec_ = GUNDAM::Nec(this->origin_data_graph_);
    this->updated_data_graph_nec_ = GUNDAM::Nec(this->updated_data_graph_);
    return;
  }

  inline bool load_data_graph_success() const {
    return this->load_data_graph_success_;
  }

 private:
  bool load_data_graph_success_;

  std::string graph_path_v_file_,
              graph_path_e_file_;

   GraphType  origin_data_graph_,
             updated_data_graph_;

  EdgeIDType origin_max_edge_id_,
            updated_max_edge_id_;

  std::unordered_map<VertexIDType,
                  std::pair<typename GUNDAM::VertexHandle<GraphType>::type,
                            typename GUNDAM::VertexHandle<GraphType>::type>>
                id_to_origin_updated_vertex_handle_;

  std::unordered_set<VertexIDType> updated_vertex_id_;

  std::vector<
  std::vector<typename GUNDAM::VertexHandle<GraphType>::type>>
                origin_data_graph_nec_, updated_data_graph_nec_;

  std::map<typename GUNDAM::VertexLabel <GraphType>::type,
           std::vector<
           typename GUNDAM::VertexHandle<GraphType>::type>> centrel_vertex_set_;
};

}; // _gar_discover

}; // grape

#endif // EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_GRAPH_PACKAGE_H_