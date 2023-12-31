#ifndef EXAMPLES_ANALYTICAL_GPARNEW_H_
#define EXAMPLES_ANALYTICAL_GPARNEW_H_

#include <string>
#include <vector>

#include "gundam/io/csvgraph.h"
//#include "gundam/graph_type/graph.h"

namespace gar {

template <class Pattern>
class GPAR {
 public:
  typedef Pattern value_type;
  using VertexType = typename Pattern::VertexType;
  using VertexIDType = typename VertexType::IDType;
  using VertexLabelType = typename VertexType::LabelType;
  using EdgeType = typename Pattern::EdgeType;
  using EdgeIDType = typename EdgeType::IDType;
  using EdgeLabelType = typename EdgeType::LabelType;
  using VertexPtr = typename Pattern::VertexPtr;
  using VertexConstPtr = typename Pattern::VertexConstPtr;
  using EdgePtr = typename Pattern::EdgePtr;
  using EdgeConstPtr = typename Pattern::EdgeConstPtr;
  using VertexSizeType = size_t;

  GPAR() = default;
  // init x_node_id = 1 and y_node_id = 2
  GPAR(const VertexLabelType x, const VertexLabelType y,
       const EdgeLabelType q) {
    this->pattern.AddVertex(1, x);
    this->pattern.AddVertex(2, y);
    this->x_node_ptr_ = this->pattern.FindConstVertex(1);
    this->y_node_ptr_ = this->pattern.FindConstVertex(2);
    this->q_edge_label_ = q;
  }
  // another constructor
  GPAR(Pattern pattern, const VertexIDType x_node_id,
       const VertexIDType y_node_id, const EdgeLabelType q) {
    this->pattern = pattern;
    this->x_node_ptr_ = this->pattern.FindConstVertex(x_node_id);
    this->y_node_ptr_ = this->pattern.FindConstVertex(y_node_id);
    this->q_edge_label_ = q;
  }
  // copy onstructor
  GPAR(const GPAR<Pattern>& b) {
    this->pattern = b.pattern;
    this->x_node_ptr_ = this->pattern.FindConstVertex(b.x_node_ptr_->id());
    this->y_node_ptr_ = this->pattern.FindConstVertex(b.y_node_ptr_->id());
    this->q_edge_label_ = b.q_edge_label_;
  }
  VertexConstPtr x_node_ptr() const { return this->x_node_ptr_; }
  VertexConstPtr y_node_ptr() const { return this->y_node_ptr_; }
  EdgeLabelType q_edge_label() const { return this->q_edge_label_; }
  bool operator<(const GPAR<Pattern>& b) const { return 1; }
  GPAR& operator=(const GPAR<Pattern>& b) {
    this->pattern = b.pattern;
    this->x_node_ptr_ = this->pattern.FindConstVertex(b.x_node_ptr_->id());
    this->y_node_ptr_ = this->pattern.FindConstVertex(b.y_node_ptr_->id());
    this->q_edge_label_ = b.q_edge_label_;
    return *this;
  }

 public:
  Pattern pattern;
  VertexConstPtr x_node_ptr_, y_node_ptr_;
  EdgeLabelType q_edge_label_;
};
template <class Pattern, class DataGraph>
class DiscoverdGPAR : public GPAR<Pattern> {
 public:
  typedef Pattern value_type;
  using VertexType = typename Pattern::VertexType;
  using VertexIDType = typename VertexType::IDType;
  using VertexLabelType = typename VertexType::LabelType;
  using EdgeType = typename Pattern::EdgeType;
  using EdgeIDType = typename EdgeType::IDType;
  using EdgeLabelType = typename EdgeType::LabelType;
  using VertexPtr = typename Pattern::VertexPtr;
  using VertexConstPtr = typename Pattern::VertexConstPtr;
  using EdgePtr = typename Pattern::EdgePtr;
  using EdgeConstPtr = typename Pattern::EdgeConstPtr;
  using VertexSizeType = size_t;
  using DataGraphVertexPtr = typename DataGraph::VertexConstPtr;
  using SuppType = std::vector<DataGraphVertexPtr>;
  using ConfType = double;

 private:
  SuppType supp_Q_;
  SuppType supp_R_;
  ConfType conf_;
  int GPAR_id;

 public:
  DiscoverdGPAR() = default;
  DiscoverdGPAR(const VertexLabelType x, const VertexLabelType y,
                const EdgeLabelType q)
      : GPAR<Pattern>(x, y, q) {
    this->conf_ = -1;
    this->supp_Q_.clear();
    this->supp_R_.clear();
    this->GPAR_id = -1;
  }
  DiscoverdGPAR(Pattern pattern, const VertexIDType x_node_id,
                const VertexIDType y_node_id, const EdgeLabelType q)
      : GPAR<Pattern>(pattern, x_node_id, y_node_id, q) {
    this->conf_ = -1;
    this->supp_Q_.clear();
    this->supp_R_.clear();
    this->GPAR_id = -1;
  }
  DiscoverdGPAR(const DiscoverdGPAR<Pattern, DataGraph>& b) : GPAR<Pattern>(b) {
    this->supp_R_ = b.supp_R_;
    this->supp_Q_ = b.supp_Q_;
    this->conf_ = b.conf_;
    this->GPAR_id = b.GPAR_id;
  }
  bool operator<(const DiscoverdGPAR<Pattern, DataGraph>& b) const { return 1; }
  DiscoverdGPAR& operator=(const DiscoverdGPAR<Pattern, DataGraph>& b) {
    this->pattern = b.pattern;
    this->x_node_ptr_ = this->pattern.FindConstVertex(b.x_node_ptr()->id());
    this->y_node_ptr_ = this->pattern.FindConstVertex(b.y_node_ptr()->id());
    this->q_edge_label_ = b.q_edge_label();
    this->supp_R_ = b.supp_R_;
    this->supp_Q_ = b.supp_Q_;
    this->conf_ = b.conf_;
    this->GPAR_id = b.GPAR_id;
    return *this;
  }
  SuppType& supp_Q() { return this->supp_Q_; }
  SuppType& supp_R() { return this->supp_R_; }
  size_t supp_Q_size() const { return this->supp_Q_.size(); }
  size_t supp_R_size() const { return this->supp_R_.size(); }
  ConfType conf() const { return this->conf_; }
  void CalConf() {
    if (!this->supp_Q_.empty()) {
      this->conf_ = static_cast<double>(this->supp_R_.size()) /
                    static_cast<double>(this->supp_Q_.size());
    }
  }
  int id() const { return this->GPAR_id; }
  void Setid(int id) { this->GPAR_id = id; }
};
template <class GPAR>
void WriteGPAR(const GPAR& gpar, uint32_t pattern_id,
               const std::string& output_dir) {
  std::string v_file =
      output_dir + "pattern_" + std::to_string(pattern_id) + "_v.csv";
  std::string e_file =
      output_dir + "pattern_" + std::to_string(pattern_id) + "_e.csv";
  GUNDAM::WriteCSVGraph(gpar.pattern, v_file, e_file);
}

}  // namespace gar

#endif