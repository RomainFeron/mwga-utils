#pragma once
#include <iomanip>
#include <sstream>
#include "cli.h"
#include "utils.h"

void output_alignability_table(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics);

void output_alignability_wig(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics);

void output_identity_wig(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics);

void output_major_allele_wig(const Parameters& parameters, const std::unordered_map<std::string, std::vector<BaseData>>& metrics);

void output_allele_count_wig(Parameters& parameters, std::unordered_map<std::string, std::vector<BaseData>> metrics);
