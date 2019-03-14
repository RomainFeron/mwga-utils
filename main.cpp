#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>


struct AlignabilityData {

    uint16_t all = 0;
    uint16_t no_N = 0;
};


int main() {

    std::string maf_file_path = "/home/romain/work/data/21genomes/roast/combined_roast.maf";
    std::ifstream maf_file(maf_file_path);

    if (not maf_file.is_open()) {
        std::cerr << "Error: input file <" << maf_file_path << "> does not exist." << std::endl;
        exit(1);
    }

    std::unordered_map<std::string, std::vector<AlignabilityData>> alignability;

    std::string line, current_field, scaffold;
    std::vector<std::string> fields(6);
    std::vector<bool> columns;
    std::vector<bool> columns_including_N;
    bool new_block = true, new_field = true;
    uint8_t field_n = 0;
    uint32_t position = 0, ref_position = 0;
    uint32_t start = 0;
    uint32_t line_count = 0;

    while(std::getline(maf_file, line)) {

        if (line.size() > 0) {

            if (line[0] == 'a') {

                new_block = true;
                scaffold = "";

            } else if (line[0] == 's') {

                field_n = 0;
                current_field = "";
                new_field = false;
                position = 0;
                ref_position = 0;

                for (auto c: line) {

                    switch (c) {

                        case ' ':

                            if (current_field != "") {
                                fields[field_n] = current_field;
                                new_field = true;
                                current_field = "";
                            }

                            break;

                        default:

                            if (new_field) {

                                ++field_n;

                                if (field_n == 6) {

                                    if (new_block) {

                                        scaffold = fields[1];
                                        start = uint(std::stoi(fields[2]));

                                        if (alignability.find(scaffold) == alignability.end()) {
                                            alignability[scaffold].reserve(uint(std::stoi(fields[5])));
                                            alignability[scaffold].resize(uint(std::stoi(fields[5])));
                                        }

                                        columns.resize(0);
                                        columns_including_N.resize(0);

                                    }
                                }

                                new_field = false;
                            }

                            if (field_n < 6) {

                                current_field += c;

                            } else {

                                if (new_block) {

                                    switch (c) {

                                        case '-':
                                            columns_including_N.push_back(false);
                                            columns.push_back(false);
                                            break;

                                        case 'N':
                                            columns_including_N.push_back(true);
                                            columns.push_back(false);
                                            break;

                                        default:
                                            columns.push_back(true);
                                            columns_including_N.push_back(true);
                                            break;
                                    }

                                } else {

                                    if (c != '-') {

                                        if (columns[position]) ++alignability[scaffold][ref_position + start].no_N;

                                        if (columns_including_N[position]) {
                                            ++alignability[scaffold][ref_position + start].all;
                                            ++ref_position;
                                        }

                                    }
                                }
                                ++position;
                            }

                            break;
                    }

                }

                if (new_block) {
                    new_block = false;
                }

            }
        }

        if (line_count % 10000 == 0 and line_count != 0) std::cout << std::setprecision(5) << "Processed <" << double(line_count) / 1000000.0 << "> M. lines" << std::endl;

        ++line_count;
    }

    std::string output_file_path = "/home/romain/work/code/wga_stats/test_cpp.tsv";
    std::ofstream output_file(output_file_path);

    if (not output_file.is_open()) {
        std::cerr << "Error: output file <" << output_file_path << "> cannot be created." << std::endl;
        exit(1);
    }

    output_file << "Scaffold\tSize\tAlignability\tBp\tBp_N\n";

    for (auto& scaffold: alignability) {

        std::vector<uint32_t> distribution(21);
        std::vector<uint32_t> distribution_including_N(21);

        for (auto nuc: scaffold.second) {
            ++distribution[nuc.all];
            ++distribution_including_N[nuc.no_N];
        }

        for (uint i=0; i<21; ++i) {
            output_file << scaffold.first << "\t" << scaffold.second.size() << "\t" << i + 1 << "\t" << distribution[i] << "\t" << distribution_including_N[i] << "\n";
        }

        distribution.resize(0);
        distribution_including_N.resize(0);

    }

    return 0;
}
