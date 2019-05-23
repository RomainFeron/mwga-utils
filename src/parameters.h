#pragma once
#include <string>
#include <vector>


// Data structure to store alignability data for a single nucleotide
struct BaseData {
    uint16_t alignability = 0;  // Number of assemblies in which this base was aligned
    float identity = 0.0;
    std::vector<uint8_t> alleles {0, 0, 0, 0, 0, 0};
    bool is_N = false;  // Boolean indicating whether this base is an N in the reference assembly
    char ref_base = ' ';
};


struct Parameters {

    // Subcommand to execute from CLI
    std::string command = "";

    // I/O parameters
    std::string maf_file_path = "";
    std::string genome_file_path = "";
    std::string alignability_table_file_path = "";
    std::string alignability_wig_file_path = "";
    std::string identity_wig_file_path = "";
    std::string major_allele_wig_file_path = "";
    std::string allele_count_wig_file_path = "";

};
