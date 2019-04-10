#include "maf_stats.h"

MafStats::MafStats(int argc, char* argv[]) {

    this->parameters = parse_args(argc, argv);

    // Open input file
    this->maf_file.open(parameters.maf_file_path);

    if (not maf_file.is_open()) {
        std::cerr << "Error opening input file <" << parameters.maf_file_path << ">." << std::endl;
        exit(1);
    }
}


void MafStats::run() {

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
                                        if (this->metrics.find(scaffold) == this->metrics.end()) {
                                            this->metrics[scaffold].reserve(uint(std::stoi(fields[5])));
                                            this->metrics[scaffold].resize(uint(std::stoi(fields[5])));
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

                                    this->metrics[scaffold][ref_position + start].ref_base = c;

                                    switch (c) {

                                        case '-':  // Gaps are excluded in all cases
                                            columns.push_back(false);
                                            break;

                                        case 'N':  // When N, update the is_N field for this base
                                            columns.push_back(true);
                                            this->metrics[scaffold][ref_position + start].is_N = true;  // ref_position is incremented for all bases except gaps
                                            ++ref_position;
                                            break;

                                        default:  // All other bases are included
                                            columns.push_back(true);
                                            ++ref_position;
                                            break;
                                    }

                                } else {

                                    if (columns[position]) {
                                        if (c != '-') {  // Gaps are excluded from non-reference assemblies too
                                            ++this->metrics[scaffold][ref_position + start].alignability;  // Increment alignability for this position
                                            this->metrics[scaffold][ref_position + start].identity += (c == this->metrics[scaffold][ref_position + start].ref_base);  // Increment identity for this position
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

    if (parameters.alignability_table_file_path != "") this->output_alignability_table();
    if (parameters.alignability_wig_file_path != "") this->output_alignability_wig();
    if (parameters.identity_wig_file_path != "") this->output_identity_wig();
}



void MafStats::output_alignability_table() {

    std::cerr << "Generating alignability distribution table..." << std::endl;

    // Open output file
    std::ofstream output_file;
    output_file.open(parameters.alignability_table_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening output file <" << parameters.alignability_table_file_path << ">." << std::endl;
        exit(1);
    }
    output_file << this->standard_output_header;

    // Output the data
    for (auto& scaffold: this->metrics) {

        std::vector<uint32_t> distribution(21);  // Distribution of alignability including Ns for the current scaffold
        std::vector<uint32_t> distribution_excluding_N(21);  // Distribution of alignability excluding Ns for the current scaffold
        uint32_t scaffold_size_excluding_N = 0;  // Total scaffold size after excluding N

        for (auto nuc: scaffold.second) {
            ++distribution[nuc.alignability];
            if (not nuc.is_N) {
                ++distribution_excluding_N[nuc.alignability];
                ++scaffold_size_excluding_N;
            }
        }

        for (uint i=0; i<21; ++i) {
            std::ostringstream output_line;
            output_line << scaffold.first << "\t" << scaffold.second.size() << "\t" << scaffold_size_excluding_N
                        << "\t" << i + 1 << "\t" << distribution[i] << "\t" << distribution_excluding_N[i] << "\n";
            output_file << output_line.str();
        }
    }
}



void MafStats::output_alignability_wig() {

    std::cerr << "Generating alignability wig file ..." << std::endl;

    // Open wig file
    std::ofstream output_file;
    output_file.open(parameters.alignability_wig_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening wig file <" << parameters.alignability_wig_file_path << ">." << std::endl;
        exit(1);
    }

    // Output the data
    for (auto& scaffold: this->metrics) {
        output_file << "fixedStep chrom=" << split(scaffold.first, ".")[1] << " start=1 step=1\n";
        for (auto nuc: scaffold.second) {
            output_file << nuc.alignability + 1 << "\n";
        }
    }
}



void MafStats::output_identity_wig() {

    std::cerr << "Generating identity wig file ..." << std::endl;

    // Open wig file
    std::ofstream output_file;
    output_file.open(parameters.identity_wig_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening wig file <" << parameters.identity_wig_file_path << ">." << std::endl;
        exit(1);
    }

    // Output the data
    for (auto& scaffold: this->metrics) {
        output_file << "fixedStep chrom=" << split(scaffold.first, ".")[1] << " start=1 step=1\n";
        for (auto nuc: scaffold.second) {
            output_file << std::setprecision(3) << (nuc.identity + 1) / (nuc.alignability + 1) << "\n";
        }
    }
}
