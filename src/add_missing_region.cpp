#include "add_missing_region.h"

void add_missing_regions(Parameters& parameters) {

    // Open input file
    std::ifstream maf_file;
    maf_file.open(parameters.maf_file_path);

    if (not maf_file.is_open()) {
        std::cerr << "Error opening input file <" << parameters.maf_file_path << ">." << std::endl;
        exit(1);
    }

    std::string line, current_field, scaffold;
    std::vector<std::string> fields(6);
    bool new_block = true, new_field = true;
    uint8_t field_n = 0;
    uint32_t start = 0, end = 0, size = 0;
    uint32_t line_count = 0;

    std::unordered_map<std::string, std::vector<uint>> coverage;

    std::unordered_map<std::string, uint> scaffold_lengths;

    std::string o_line = "";

    std::cerr << "Processing MAF file ..." << std::endl;

    // Read the MAF file
    while(std::getline(maf_file, line)) {

        o_line = "";

        if (line.size() > 0) {

            if (line[0] == 'a') {  // Lines starting with 'a' indicate the start of a new block

                new_block = true;  // Record that a new block has started
                scaffold = "";
                o_line = line;

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

                            o_line += c;

                            break;

                        default:

                            if (field_n < 6) {

                                // Fields : 's', scaffold, start_pos, alignment_length, strand, scaffold_length, sequence
                                // The sequence is not stored in current field, it is processed separately

                                current_field += c;  // Update current field string except for the sequence field

                                if (field_n < 5) o_line += c; else o_line += static_cast<char>(std::toupper(c));

                            } else {
                                o_line += static_cast<char>(std::toupper(c));
                            }

                            if (new_field) { // Check if this is the start of a field

                                ++field_n;

                                if (field_n == 6 and new_block) {  // Check if this is the first sequence in a block

                                    scaffold = fields[1];
                                    start = uint(std::stoi(fields[2]));
                                    end = start + uint(std::stoi(fields[3]));

                                    // Initialize a vector of the right size for the scaffold if this is the first time this scaffold is encountered
                                    if (coverage.find(scaffold) == coverage.end()) {
                                        size = static_cast<uint32_t>(std::stoi(fields[5]));
                                        coverage[scaffold] = std::vector<uint>(size, 0);
                                        scaffold_lengths[scaffold] = size;
                                    }

                                    for (auto i=start; i<end; ++i) ++coverage[scaffold][i];

                                }

                                new_field = false;
                            }

                            break;
                    }
                }

                if (new_block) {
                    new_block = false;
                }

            } else {
                o_line = line;
            }
        }

        if (line_count % 1000000 == 0 and line_count != 0) std::cerr << "  - Processed " << line_count << " lines" << std::endl;

        ++line_count;

        std::cout << o_line << "\n";
    }


    std::unordered_map<std::string, std::string> genome;
    std::ifstream genome_file;
    genome_file.open(parameters.genome_file_path);

    std::string fa_header = "", fa_sequence = "";
    std::vector<std::string> tmp;
    while (std::getline(genome_file, line)) {
        if (line[0] == '>') {
            if (fa_header != "") genome[fa_header] = fa_sequence;
            fa_header = line.substr(1);
            tmp = split(fa_header, ":");
            fa_header = tmp[0] + "." + tmp[1];
            fa_sequence = "";
        } else {
            fa_sequence += line;
        }
    }
    genome[fa_header] = fa_sequence;

    for (auto s: genome) std::cerr << s.first << " : " << s.second.size() << std::endl;


    std::string sequence = "";
    uint seq_start = 0, seq_length = 0;
    bool new_sequence = true;
    for (auto& scaffold: coverage) {
        sequence = "";
        seq_start = 0;
        seq_length = 0;
        new_sequence = true;
        for (uint i=0; i<scaffold.second.size(); ++i) {
            if (scaffold.second[i] == 1) {
                if (not new_sequence) {
                    sequence = genome[scaffold.first].substr(seq_start, seq_length);
                    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                    std::cout << "a score=NA\ns " << scaffold.first << " " << seq_start << " " << seq_length << " + " << scaffold_lengths[scaffold.first] << " " << sequence << "\n\n";
                }
                sequence = "";
                seq_start = 0;
                seq_length = 0;
                new_sequence = true;
            } else {
                if (new_sequence) {
                    seq_start = i;
                    new_sequence = false;
                }
                ++seq_length;
            }
        }

        if (not new_sequence) {
            sequence = genome[scaffold.first].substr(seq_start, seq_length);
            std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
            std::cout << "a score=NA\ns " << scaffold.first << " " << seq_start << " " << seq_length << " + " << scaffold_lengths[scaffold.first] << " " << sequence << "\n\n";
        }

    }

    maf_file.close();
}
