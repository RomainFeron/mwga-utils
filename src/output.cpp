#include "output.h"


void output_alignability_table(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics) {

    std::cerr << "Generating alignability distribution table..." << std::endl;

    // Open output file
    std::ofstream output_file;
    output_file.open(parameters.alignability_table_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening output file <" << parameters.alignability_table_file_path << ">." << std::endl;
        exit(1);
    }

    // Output the header
    output_file << "Scaffold\tSize\tSize_no_N\tAlignability\tCount\tCount_no_N\n";

    std::ostringstream output_line;

    // Output the data
    for (auto& scaffold: metrics) {

        std::vector<uint32_t> distribution(22);  // Distribution of alignability including Ns for the current scaffold
        std::vector<uint32_t> distribution_excluding_N(22);  // Distribution of alignability excluding Ns for the current scaffold
        uint32_t scaffold_size_excluding_N = 0;  // Total scaffold size after excluding N

        for (auto& nuc: scaffold.second) {
            ++distribution[nuc.alignability];
            if (not nuc.is_N) {
                ++distribution_excluding_N[nuc.alignability];
                ++scaffold_size_excluding_N;
            }
        }

        for (uint i=0; i<22; ++i) {
            output_line.str("");
            output_line.clear();
            output_line << scaffold.first << "\t" << scaffold.second.size() << "\t" << scaffold_size_excluding_N
                        << "\t" << i + 1 << "\t" << distribution[i] << "\t" << distribution_excluding_N[i];
            output_file << output_line.str() << "\n";
        }

    }

    output_file.close();
}



void output_alignability_wig(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics) {

    std::cerr << "Generating alignability wig file ..." << std::endl;

    // Open wig file
    std::ofstream output_file;
    output_file.open(parameters.alignability_wig_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening wig file <" << parameters.alignability_wig_file_path << ">." << std::endl;
        exit(1);
    }

    // Output the data
    for (auto& scaffold: metrics) {
        output_file << "fixedStep chrom=" << split(scaffold.first, ".")[1] << " start=1 step=1\n";
        for (auto& nuc: scaffold.second) {
            output_file << nuc.alignability + 1 << "\n";
        }
    }
}



void output_identity_wig(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics) {

    std::cerr << "Generating identity wig file ..." << std::endl;

    // Open wig file
    std::ofstream output_file;
    output_file.open(parameters.identity_wig_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening wig file <" << parameters.identity_wig_file_path << ">." << std::endl;
        exit(1);
    }

    // Output the data
    for (auto& scaffold: metrics) {
        output_file << "fixedStep chrom=" << split(scaffold.first, ".")[1] << " start=1 step=1\n";
        for (auto& nuc: scaffold.second) {
            output_file << std::setprecision(3) << (nuc.identity + 1) / (nuc.alignability + 1) << "\n";
        }
    }
}



void output_major_allele_wig(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics) {

    std::cerr << "Generating major allele wig file ..." << std::endl;

    // Open wig file
    std::ofstream output_file;
    output_file.open(parameters.major_allele_wig_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening wig file <" << parameters.major_allele_wig_file_path << ">." << std::endl;
        exit(1);
    }

    float major_allele = 0.0, total_count = 0.0;

    // Output the data
    for (auto& scaffold: metrics) {
        output_file << "fixedStep chrom=" << split(scaffold.first, ".")[1] << " start=1 step=1\n";
        for (auto& nuc: scaffold.second) {
            total_count = static_cast<float>(nuc.alleles[0]) + static_cast<float>(nuc.alleles[1]) + static_cast<float>(nuc.alleles[2]) +
                          static_cast<float>(nuc.alleles[3]) + static_cast<float>(nuc.alleles[4]) + static_cast<float>(nuc.alleles[5]);
            major_allele = static_cast<float>(*std::max_element(nuc.alleles.begin(), nuc.alleles.end()));
            output_file << std::setprecision(3) << major_allele / total_count << "\n";
        }
    }
}



void output_allele_count_wig(Parameters& parameters, std::unordered_map<std::string, std::vector<BaseData>> metrics) {

    std::cerr << "Generating allele count wig file ..." << std::endl;

    // Open wig file
    std::ofstream output_file;
    output_file.open(parameters.allele_count_wig_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening wig file <" << parameters.allele_count_wig_file_path << ">." << std::endl;
        exit(1);
    }

    uint allele_count = 0;

    // Output the data
    for (auto& scaffold: metrics) {
        output_file << "fixedStep chrom=" << split(scaffold.first, ".")[1] << " start=1 step=1\n";
        for (auto& nuc: scaffold.second) {

            switch (nuc.ref_base) {

                case 'A':
                    ++nuc.alleles[0];
                    break;

                case 'T':
                    ++nuc.alleles[1];
                    break;

                case 'G':
                    ++nuc.alleles[2];
                    break;

                case 'C':
                    ++nuc.alleles[3];
                    break;

                case 'N':
                    ++nuc.alleles[4];
                    break;

                case '-':
                    ++nuc.alleles[5];
                    break;
            }

            allele_count = (nuc.alleles[0] > 0) + (nuc.alleles[1] > 0) + (nuc.alleles[2] > 0) + (nuc.alleles[3] > 0) + (nuc.alleles[4] > 0) + (nuc.alleles[5] > 0);
            output_file << allele_count << "\n";
        }
    }
}
