#pragma once
#include <iomanip>
#include <iostream>
#include <thread>
#include "maf_parser.h"
#include "utils.h"


typedef std::unordered_map<std::string, uint> BpAlignedData;

void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, BpAlignedData& ref_coverage, bool& parsing_ended, unsigned long batch_size);
std::ofstream open_ofile(std::string& output_file_path);
