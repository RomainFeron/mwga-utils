#pragma once
#include <string>

struct Parameters {

    // I/O parameters
    std::string maf_file_path = "";
    std::string alignability_table_file_path = "";
    std::string alignability_wig_file_path = "";
    std::string identity_wig_file_path = "";
    std::string major_allele_wig_file_path = "";
    std::string allele_count_wig_file_path = "";

};
