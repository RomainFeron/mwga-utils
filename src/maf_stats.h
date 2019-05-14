#pragma once
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "parser.h"
#include "utils.h"


// Data structure to store alignability data for a single nucleotide
struct BaseData {
    uint16_t alignability = 0;  // Number of assemblies in which this base was aligned
    float identity = 0.0;
    uint8_t alleles[6] {0, 0, 0, 0, 0, 0};
    bool is_N = false;  // Boolean indicating whether this base is an N in the reference assembly
    char ref_base = ' ';
};


class MafStats {

    public:

        Parameters parameters;
        std::ifstream maf_file;
        std::ofstream output_file;
        std::ofstream wig_output_file;

        std::string standard_output_header = "Scaffold\tSize\tSize_no_N\tAlignability\tCount\tCount_no_N\n";

        std::unordered_map<std::string, std::vector<BaseData>> metrics;  // Store alignability results

        MafStats();
        MafStats(int argc, char* argv[]);
        void run();
        void output_alignability_table();
        void output_alignability_wig();
        void output_identity_wig();
        void output_major_allele_wig();
        void output_allele_count_wig();
};

