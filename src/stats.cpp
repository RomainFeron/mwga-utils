#include "bp_aligned.h"

void bp_aligned(Parameters& parameters) {

    // Open input file
    std::ifstream maf_file;
    maf_file.open(parameters.maf_file_path);

    if (not maf_file.is_open()) {
        std::cerr << "Error opening input file <" << parameters.maf_file_path << ">." << std::endl;
        exit(1);
    }

    std::string line, current_field, species;
    std::vector<std::string> fields(6);
    bool new_block = true, new_field = true;
    uint8_t field_n = 0;
    uint32_t block_size = 0, size = 0;
    uint32_t line_count = 0;

    std::unordered_map<std::string, uint> coverage;

    std::cerr << "Processing MAF file ..." << std::endl;

    // Read the MAF file
    while(std::getline(maf_file, line)) {

        if (line.size() > 0) {

            if (line[0] == 'a') {  // Lines starting with 'a' indicate the start of a new block

                new_block = true;  // Record that a new block has started

            } else if (line[0] == 's') {  // Lines starting with 's' are sequence lines inside a block

                field_n = 0;
                current_field = "";
                new_field = false;

                for (auto c: line) {  // Iterate over the line

                    switch (c) {

                        case ' ':  // Fields are delimited by varying number of spaces

                            if (current_field != "") {  // Checking that the previous character was not a space (i.e this is the first space in the delimiter).
                                fields[field_n] = current_field;
                                new_field = true;  // Record that a new field has started
                                current_field = "";  // Current field is reset to "" when a field ends
                            }

                            break;

                        default:

                            if (field_n < 6) {

                                // Fields : 's', scaffold, start_pos, alignment_length, strand, scaffold_length, sequence
                                // The sequence is not stored in current field, it is processed separately

                                current_field += c;  // Update current field string except for the sequence field

                            }

                            if (new_field) { // Check if this is the start of a field

                                ++field_n;

                                if (field_n == 6) {  // Check if this is the first sequence in a block

                                    species = split(fields[1], ".")[0];
                                    block_size = uint(std::stoi(fields[3]));
                                    coverage[species] += block_size;

                                }

                                new_field = false;
                            }

                            break;
                    }
                }

                if (new_block) {
                    new_block = false;
                }

            }
        }

        if (line_count % 1000000 == 0 and line_count != 0) std::cerr << "  - Processed " << line_count << " lines" << std::endl;

        ++line_count;
    }

    for (auto& species: coverage) {
        std::cout << species.first << "\t" << species.second << std::endl;
    }

    maf_file.close();
}
