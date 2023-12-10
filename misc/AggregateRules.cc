#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>

using namespace std;

int main (int argc, char** argv) {
  string input_file_prefix(argv[1]);
  string input_level_num_str(argv[2]);
  string input_worker_num_str(argv[3]);
  string output_file_prefix(argv[4]);

  int input_level_num = stoi(input_level_num_str);
  int input_worker_num = stoi(input_worker_num_str);

  int global_gar_idx = 0;

  string output_t_file_str  = output_file_prefix + "/gar_t_set.csv";
  string output_v_file_str  = output_file_prefix + "/gar_v_set.csv";
  string output_e_file_str  = output_file_prefix + "/gar_e_set.csv";
  string output_x_file_str  = output_file_prefix + "/gar_x_set.csv";
  string output_y_file_str  = output_file_prefix + "/gar_y_set.csv";

  ofstream output_t_file(output_t_file_str), output_v_file(output_v_file_str),
           output_e_file(output_e_file_str), output_x_file(output_x_file_str),
           output_y_file(output_y_file_str);


  for (int level_num = 2; level_num < input_level_num; level_num++) {
    for (int worker_id = 0; worker_id < input_worker_num; worker_id++) {
      std::cout << "here" << level_num << " " << worker_id << std::endl;
      string current_t_file_str = input_file_prefix + "/gar_level_" + std::to_string(level_num)
                                    + "_worker_" + std::to_string(worker_id)
                                    + "_t_set.csv";
      string current_v_file_str = input_file_prefix + "/gar_level_" + std::to_string(level_num)
                                    + "_worker_" + std::to_string(worker_id)
                                    + "_v_set.csv";
      string current_e_file_str = input_file_prefix + "/gar_level_" + std::to_string(level_num)
                                    + "_worker_" + std::to_string(worker_id)
                                    + "_e_set.csv";
      string current_x_file_str = input_file_prefix + "/gar_level_" + std::to_string(level_num)
                                    + "_worker_" + std::to_string(worker_id)
                                    + "_x_set.csv";
      string current_y_file_str = input_file_prefix + "/gar_level_" + std::to_string(level_num)
                                    + "_worker_" + std::to_string(worker_id)
                                    + "_y_set.csv";
      ifstream input_t_file(current_t_file_str), input_v_file(current_v_file_str),
               input_e_file(current_e_file_str), input_x_file(current_x_file_str),
               input_y_file(current_y_file_str);

      if (!input_t_file.is_open() || !input_v_file.is_open() || !input_e_file.is_open() || !input_x_file.is_open() || !input_y_file.is_open()) {
        std::cout << "can not open input t file " << current_t_file_str << std::endl;
        return -1;
      }
      std::unordered_map<int, int> local_id_to_global_id;
      local_id_to_global_id[-1] = -1;
      int t_worker_id, t_v_num, t_gar_id, t_parent_id, t_supp;
      double probability;
      std::string str, match_file_name;
      while (getline(input_t_file, str)) {
        stringstream ss(str);
        ss >> t_worker_id >> t_v_num >> t_gar_id >> t_parent_id >> t_supp;
        ss >> match_file_name;
        ss >> probability;

        local_id_to_global_id[t_gar_id] = global_gar_idx;
        output_t_file << t_worker_id << " " << t_v_num 
                      << " " << global_gar_idx << " " << local_id_to_global_id[t_parent_id]
                      << " " << t_supp << " " << match_file_name 
                      << " " << probability << "\n";
        global_gar_idx++;
      }
      input_t_file.close();

      getline(input_v_file, str);
      if (level_num == 2 && worker_id == 0) {
        output_v_file << str << "\n";
      }
      while (getline(input_v_file, str)) {
        std::stringstream ss(str);
        string tmp_str;
        unsigned index = 0;
        while (getline(ss, tmp_str, ',')) {
          if (index != 2) {
            output_v_file << tmp_str << ",";
            index++;
          } else {
            int local_gar_idx = stoi(tmp_str);
            int global_gar_idx = local_id_to_global_id[local_gar_idx];
            output_v_file << std::to_string(global_gar_idx) << "\n";
            index++;
          }
          if (index > 3) {
            std::cout << "error v index " << std::endl;
            return -1;
          }
        }
      }
      input_v_file.close();

      getline(input_e_file, str);
      if (level_num == 2 && worker_id == 0) {
        output_e_file << str << "\n";
      }
      while (getline(input_e_file, str)) {
        std::stringstream ss(str);
        string tmp_str;
        unsigned index = 0;
        while (getline(ss, tmp_str, ',')) {
          if (index != 4) {
            output_e_file << tmp_str << ",";
            index++;
          } else {
            int local_gar_idx = stoi(tmp_str);
            int global_gar_idx = local_id_to_global_id[local_gar_idx];
            output_e_file << std::to_string(global_gar_idx) << "\n";
            index++;
          }
          if (index > 5) {
            std::cout << "error e index " << std::endl;
            return -1;
          }
        }
      }
      input_e_file.close();

      getline(input_x_file, str);
      if (level_num == 2 && worker_id == 0) {
        output_x_file << str << "\n";
      }
      while (getline(input_x_file,str)) {
        std::stringstream ss(str);
        string tmp_str;
        unsigned index = 0;
        while (getline(ss, tmp_str, ',')) {
          if (index != 7) {
            output_x_file << tmp_str << ",";
            index++;
          } else {
            int local_gar_idx = stoi(tmp_str);
            int global_gar_idx = local_id_to_global_id[local_gar_idx];
            output_x_file << std::to_string(global_gar_idx) << "\n";
            index++;
          }
        }
      }
      input_x_file.close();

      getline(input_y_file, str);
      if (level_num == 2 && worker_id == 0) {
        output_y_file << str << "\n";
      }
      while (getline(input_y_file, str)) {
        std::stringstream ss(str);
        string tmp_str;
        unsigned index = 0;
        while (getline(ss, tmp_str, ',')) {
          if (index != 7) {
            output_y_file << tmp_str << ",";
            index++;
          } else {
            int local_gar_idx = stoi(tmp_str);
            int global_gar_idx = local_id_to_global_id[local_gar_idx];
            output_y_file << std::to_string(global_gar_idx) << "\n";
            index++;
          }
        }
      }
      input_y_file.close();
    }
  }
  output_t_file.close();
  output_v_file.close();
  output_x_file.close();
  output_y_file.close();
  output_e_file.close();
  return 0;
}
