#ifndef EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_SERIALIZE_H_
#define EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_SERIALIZE_H_

#include "util/log.h"

namespace grape {

namespace _gar_discover{

template <typename DataGraphLiteralType>
class UpdateLiteral {
 public:
  inline const int id() const{
    return this->id_;
  }

  inline void set_id(int id) {
    this->id_ = id;
    return;
  }

  inline std::vector<DataGraphLiteralType>& rhs_literals() {
    return this->rhs_literals_;
  }

  inline const std::vector<DataGraphLiteralType>& rhs_literals() const {
    return this->rhs_literals_;
  }

  inline std::vector<DataGraphLiteralType>& lhs_literals() {
    return this->lhs_literals_;
  }

  inline const std::vector<DataGraphLiteralType>& lhs_literals() const {
    return this->lhs_literals_;
  }

 private:
  int id_;

  std::vector<DataGraphLiteralType> rhs_literals_,
                                    lhs_literals_;
};

template <typename DataGraphLiteralType>
std::string& operator<<(
    std::string& out_string,
    const UpdateLiteral<DataGraphLiteralType>& update_literal) {
  out_string = std::move(out_string) + " <update_literal";

  out_string = std::move(out_string)
             + " " + GUNDAM::ToString(update_literal.id());

  // serilize possible rhs_literals set
  out_string = std::move(out_string) + " <rhs_literals";
  for (auto literal_cit  = update_literal.rhs_literals().begin();
            literal_cit != update_literal.rhs_literals(). end ();
            literal_cit++){
    out_string = std::move(out_string) + " ";
    out_string << *literal_cit;
  }
  out_string = std::move(out_string) + " >"; // rhs_literals

  // serilize all lhs_literals set
  out_string = std::move(out_string) + " <lhs_literals";
  for (auto literal_cit  = update_literal.lhs_literals().begin();
            literal_cit != update_literal.lhs_literals(). end ();
            literal_cit++){
    out_string = std::move(out_string) + " ";
    out_string << *literal_cit;
  }
  out_string = std::move(out_string) + " >"; // lhs_literals

  out_string = std::move(out_string) + " >";  // update_literal

  return out_string;
}

template <typename DataGraphLiteralType>
std::string& operator>>(
    std::string& in_string,
    UpdateLiteral<DataGraphLiteralType>& update_literal) {

  std::stringstream ss;
  ss << in_string;

  std::string str;

  ss >> str;
  assert(str == "<update_literal");

  // deserilize expand node id
  ss >> str;
  int expand_node_id = GUNDAM::StringToDataType<int>(str);

  update_literal.set_id(expand_node_id);

  update_literal.rhs_literals().clear();
  // deserilize possible literal set
  ss >> str;
  assert(str == "<rhs_literals");

  while (ss.peek() == ' ') {
    char c = ss.get();
  }

  while (ss.peek() != '>') {
    DataGraphLiteralType literal;
    getline(ss, str);
    str >> literal;
    ss.clear();
    ss << str;
    update_literal.rhs_literals().emplace_back(literal);
    while (ss.peek() == ' ') {
      char c = ss.get();
    }
  }
  ss >> str;
  assert(str == ">");  // rhs_literals

  update_literal.lhs_literals().clear();
  // deserilize possible literal set
  ss >> str;
  assert(str == "<lhs_literals");

  while (ss.peek() == ' ') {
    char c = ss.get();
  }

  while (ss.peek() != '>') {
    DataGraphLiteralType literal;
    getline(ss, str);
    str >> literal;
    ss.clear();
    ss << str;
    update_literal.lhs_literals().emplace_back(literal);
    while (ss.peek() == ' ') {
      char c = ss.get();
    }
  }

  ss >> str;
  assert(str == ">");  // lhs_literals

  ss >> str;
  assert(str == ">");  // update_literal

  getline(ss, in_string);
  if (ss.fail()) 
    in_string.clear();
    
  return in_string;
}

template <typename GraphPatternType>
std::string& operator<<(
    std::string& out_string,
    const std::tuple<GraphPatternType, size_t, bool>& deliver_pattern) {
  out_string = std::move(out_string) + " <pattern";

  std::string str;

  str << std::get<0>(deliver_pattern);

  out_string = std::move(out_string) + str;

  out_string = std::move(out_string) + " " + std::to_string(std::get<1>(deliver_pattern));

  out_string = std::move(out_string) + " " + std::to_string(std::get<2>(deliver_pattern));

  out_string = std::move(out_string) + " >"; // pattern

  return out_string;
}

template <typename GraphPatternType>
std::string& operator>>(
    std::string& in_string,
    std::tuple<GraphPatternType, size_t, bool>& deliver_pattern) {
  using namespace GUNDAM;

  std::stringstream ss;
  ss << in_string;

  std::string str;

  ss >> str;
  assert(str == "<pattern");

  // deserilize graph pattern
  getline(ss, str);
  str >> std::get<0>(deliver_pattern);
  ss.clear();
  ss << str;

  // pattern id 
  ss >> str;
  std::get<1>(deliver_pattern) = std::stoi(str);

  // pattern only verify is legal
  ss >> str;
  std::get<2>(deliver_pattern) = std::stoi(str);

  ss >> str;
  assert(str == ">");  // pattern

  getline(ss, in_string);
  if (ss.fail()) 
    in_string.clear();

  return in_string;
}

}; // namespace _gar_discover

}  // namespace grape

#endif // EXAMPLES_ANALYTICAL_APPS_GAR_DISCOVER_SERIALIZE_H_