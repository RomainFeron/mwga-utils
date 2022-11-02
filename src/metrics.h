#pragma once
#include <iomanip>
#include <iostream>
#include <thread>
#include "maf_parser.h"
#include "utils.h"

struct Position {
    uint alignability = 0;
    uint identity = 0;
};

typedef std::unordered_map<std::string, std::vector<Position>> Metrics;

void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, std::mutex& results_mutex, bool& parsing_ended, Metrics& metrics, unsigned long batch_size = 100);

std::ofstream open_ofile(std::string& output_file_path);
