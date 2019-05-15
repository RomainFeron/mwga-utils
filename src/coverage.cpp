#include "coverage.h"

void check_coverage(Parameters& parameters) {

    // Open input file
    std::ifstream maf_file;
    maf_file.open(parameters.maf_file_path);

    if (not maf_file.is_open()) {
        std::cerr << "Error opening input file <" << parameters.maf_file_path << ">." << std::endl;
        exit(1);
    }

    std::string line, current_field, scaffold;
    std::vector<std::string> fields(6);
    std::vector<bool> columns;
    bool new_block = true, new_field = true;
    uint8_t field_n = 0;
    uint32_t position = 0, ref_position = 0;
    uint32_t start = 0;
    uint32_t line_count = 0;

    std::cerr << "Processing MAF file ..." << std::endl;

    // Read the MAF file
    while(std::getline(maf_file, line)) {

        if (line.size() > 0) {

            if (line[0] == 'a') {  // Lines starting with 'a' indicate the start of a new block

                new_block = true;  // Record that a new block has started
                scaffold = "";

            } else if (line[0] == 's') {  // Lines starting with 's' are sequence lines inside a block

                field_n = 0;
                current_field = "";
                new_field = false;
                position = 0;
                ref_position = 0;

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

                            if (new_field) { // Check if this is the start of a field

                                ++field_n;

                                if (field_n == 6) {

                                    if (new_block) {  // Check if this is the first sequence in a block

                                        scaffold = fields[1];
                                        start = uint(std::stoi(fields[2]));

                                        // Initialize a vector of the right size for the scaffold if this is the first time this scaffold is encountered
                                        if (metrics.find(scaffold) == metrics.end()) {
                                            metrics[scaffold].reserve(uint(std::stoi(fields[5])));
                                            metrics[scaffold].resize(uint(std::stoi(fields[5])));
                                        }

                                        // Reset the columns vector that store indices for non-gap positions in the reference assembly
                                        columns.resize(0);  // All positions excluding gaps

                                    }
                                }

                                new_field = false;
                            }

                            if (field_n < 6) {

                                // Fields : 's', scaffold, start_pos, alignment_length, strand, scaffold_length, sequence
                                // The sequence is not stored in current field, it is processed separately

                                current_field += c;  // Update current field string except for the sequence field

                            } else {

                                if (new_block) {

                                    metrics[scaffold][ref_position + start].ref_base = c;

                                    switch (c) {

                                        case 'A':
                                            ++metrics[scaffold][ref_position + start].alleles[0];
                                            columns.push_back(true);
                                            ++ref_position;
                                            break;

                                        case 'T':
                                            ++metrics[scaffold][ref_position + start].alleles[1];
                                            columns.push_back(true);
                                            ++ref_position;
                                            break;

                                        case 'G':
                                            ++metrics[scaffold][ref_position + start].alleles[2];
                                            columns.push_back(true);
                                            ++ref_position;
                                            break;

                                        case 'C':
                                            ++metrics[scaffold][ref_position + start].alleles[3];
                                            columns.push_back(true);
                                            ++ref_position;
                                            break;

                                        case 'N':  // When N, update the is_N field for this base
                                            columns.push_back(true);
                                            metrics[scaffold][ref_position + start].is_N = true;  // ref_position is incremented for all bases except gaps
                                            ++metrics[scaffold][ref_position + start].alleles[4];
                                            ++ref_position;
                                            break;

                                        case '-':  // Gaps are excluded in all cases
                                            columns.push_back(false);
                                            break;

                                        default:
                                            break;
                                    }

                                } else {

                                    if (columns[position]) {

                                        if (c != '-') {  // Gaps are excluded from non-reference assemblies too
                                            ++metrics[scaffold][ref_position + start].alignability;  // Increment alignability for this position
                                            metrics[scaffold][ref_position + start].identity += (c == metrics[scaffold][ref_position + start].ref_base);  // Increment identity for this position
                                        }

                                        switch (c) {

                                            case 'A':
                                                ++metrics[scaffold][ref_position + start].alleles[0];
                                                break;

                                            case 'T':
                                                ++metrics[scaffold][ref_position + start].alleles[1];
                                                break;

                                            case 'G':
                                                ++metrics[scaffold][ref_position + start].alleles[2];
                                                break;

                                            case 'C':
                                                ++metrics[scaffold][ref_position + start].alleles[3];
                                                break;

                                            case 'N':
                                                ++metrics[scaffold][ref_position + start].alleles[4];
                                                break;

                                            case '-':
                                                ++metrics[scaffold][ref_position + start].alleles[5];
                                                break;
                                        }

                                        ++ref_position;
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

        if (line_count % 1000000 == 0 and line_count != 0) std::cerr << "  - Processed " << line_count / 1000000 << " M. lines" << std::endl;

        ++line_count;
    }

    maf_file.close();
}
