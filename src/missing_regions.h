#pragma once
#include <iomanip>
#include <iostream>
#include <thread>
#include "maf_parser.h"
#include "utils.h"


typedef std::unordered_map<std::string, std::vector<bool>> MafCoverage;
typedef std::unordered_map<std::string, std::string> Assembly;

void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, MafCoverage& maf_coverage, bool& parsing_ended, ulong batch_size = 100);
