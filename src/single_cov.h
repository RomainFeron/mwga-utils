#pragma once
#include <iomanip>
#include <iostream>
#include <thread>
#include "maf_parser.h"
#include "utils.h"


typedef std::unordered_map<std::string, std::vector<uint>> RefCoverage;

void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, RefCoverage& ref_coverage, bool& parsing_ended, unsigned long batch_size);
