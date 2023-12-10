#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace std;

int main (int argc, char** argv) {
  std::string input_t_file_address(argv[1]);
  std::string input_inc_e_file_address(argv[2]);

  int worker_id = 0, v_num = 0, gar_idx = 0, parent_idx = 0, supp = 0;
  std::string match_file_address;
  ifstream input_t_file(input_t_file_address);
  if (!input_t_file.is_open()) {
    std::cout << "can not open input t file " << input_t_file_address << std::endl;
    return -1;
  }

  std::unordered_map<int, string> gar_id_to_file;
  std::unordered_set<unsigned> deleted_e_ids;

  double probability;
  while (input_t_file >> worker_id >> v_num >> gar_idx >> parent_idx >> supp >> match_file_address >> probability) {
    gar_id_to_file[gar_idx] = match_file_address;
  }

  ifstream input_inc_e_file(input_inc_e_file_address);
  if (!input_inc_e_file.is_open()) {
    std::cout << "can not open input_e_file " << input_inc_e_file_address << std::endl;
    return -1;
  }

  unsigned e_id, src_id, dst_id, lable, inc_state;
  std::string str;
  getline (input_inc_e_file, str);
  std::cout << "title " << str << std::endl;
  while (getline(input_inc_e_file, str)) {
    std::stringstream ss(str);
    std::string tmp_str;
    unsigned index = 0;
    while (getline(ss, tmp_str, ',')) {
      if (index == 0) {
        e_id = stoi(tmp_str);
      } else if (index == 4) {
        inc_state = stoi(tmp_str);
      }
      index++;
    }
    if (inc_state == 0) {
      deleted_e_ids.insert(e_id);
    }
  }
//  while (input_inc_e_file >> e_id >> src_id >> dst_id >> lable >> inc_state) {
//    std::cout << "inc state " << inc_state << endl;
//    if (inc_state == 0) {
//      deleted_e_ids.insert(e_id);
//    }
//  }
  std::cout << "e ids size " << deleted_e_ids.size() << std::endl;
  unsigned affected_pivots_num = 0, total_pivots_num = 0;

  for (auto &p : gar_id_to_file) {
    int gar_id = p.first;
    string match_file_name =  p.second;
//    std::cout << "match_file_name " << match_file_name << std::endl;
    ifstream match_file(match_file_name);
    if (!match_file.is_open()) {
      std::cout << "can not open input file " << std::endl;
    }
    int pivot_x_id = -1, pivot_y_id = -1;
    match_file >> pivot_x_id >> pivot_y_id;

//    std::cout << "pivot x " << pivot_x_id << " pivot y " << pivot_y_id << std::endl;

    std::vector<std::pair<unsigned, unsigned>> affected_data_pivots, unaffected_data_pivots;

    string str;

    while (getline(match_file, str)) {
      stringstream ss(str);
      unsigned data_x_id, data_y_id;
      ss >> data_x_id >> data_y_id;
      total_pivots_num++;
      bool affected = false;
      while (ss.rdbuf()->in_avail() != 0) {
        unsigned match_eid = 0;
        ss >> match_eid;
//        cout << "match eid " << match_eid << std::endl;
        if (deleted_e_ids.find(match_eid) == deleted_e_ids.end()) {
          continue;
        } else {
          affected = true;
          affected_pivots_num++;
          break;
        }
      }
      if (affected) {
        affected_data_pivots.emplace_back(std::make_pair(data_x_id, data_y_id));
      } else {
        unaffected_data_pivots.emplace_back(std::make_pair(data_x_id, data_y_id));
      }
    }
    match_file.close();
    string output_file_address = match_file_name + ".affected";
    ofstream output_file(output_file_address);
    if (!output_file.is_open()) {
      std::cout << "can not open output file " << output_file_address << std::endl;
      return -1;
    }
    for (auto &affected_p : affected_data_pivots) {
      output_file << affected_p.first << " " << affected_p.second << " 1\n";
    }
    for (auto &unaffected_p : unaffected_data_pivots) {
      output_file << unaffected_p.first << " " << unaffected_p.second << " 0\n";
    }
    output_file.close();
  }
  std::cout << "gars num " << gar_id_to_file.size() << std::endl;
  std::cout << "total pivots num " << total_pivots_num << " affected_pivots_num " << affected_pivots_num << endl;
  return 0;
}
