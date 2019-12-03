#pragma once
#include <mutex>
#include <queue>
#include "docopt/docopt.h"
#include "utils.h"


struct MafRecord {
    std::string source = "";
    uint start = 0;
    uint length = 0;
    char strand = '+';
    uint source_length = 0;
    std::string sequence = "";
};

struct MafBlock {
    std::string score = "";
    std::vector<MafRecord> records;
    uint n_records = 0;
};

typedef std::queue<MafBlock> BlocksQueue;

void maf_parser(std::ifstream& maf_file, BlocksQueue& blocks_queue, std::mutex& queue_mutex, bool& parsing_ended, uint& n_assemblies);

std::vector<MafBlock> get_batch(BlocksQueue& blocks_queue, std::mutex& queue_mutex, ulong batch_size=1000);
