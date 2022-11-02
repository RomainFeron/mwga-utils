#include "single_cov.h"

static const char USAGE[] =
R"(Check that the each sequence in the reference species is only covered once in a MAF.

    Usage:
      single_coverage <maf_file> [-t <threads>]

    Options:
      <maf_file>       Path to a MAF file.
      -h --help        Show this screen.
)";


void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, RefCoverage& ref_coverage, bool& parsing_ended, unsigned long batch_size) {

    size_t sep_pos = 0;
    std::string scaffold = "";
    std::vector<MafBlock> batch;
    bool keep_going = true;

    while (keep_going) {
        batch = get_batch(blocks_queue, queue_mutex, batch_size);
        if (batch.size() > 0) {
            for (auto block: batch) {
                sep_pos = block.records[0].source.find(".");
                scaffold = block.records[0].source.substr(sep_pos + 1);
                if (ref_coverage.find(scaffold) == ref_coverage.end()) ref_coverage[scaffold] = std::vector<uint>(block.records[0].source_length, 0);
                for (uint i=0; i<block.records[0].length; ++i) {
                    ++ref_coverage[scaffold][block.records[0].start + i];
                }
            }
        }
        queue_mutex.lock();
        if (parsing_ended and blocks_queue.size() == 0) keep_going = false;
        queue_mutex.unlock();
    }
}



int main(int argc, char *argv[]) {

    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "Parser 0.1");
    std::string maf_file_path = args["<maf_file>"].asString();

    RefCoverage ref_coverage;

    // Process MAF file
    std::cerr << "Processing the MAF file" << std::endl;
    std::ifstream maf_file = check_open(maf_file_path);
    bool parsing_ended = false;

    BlocksQueue blocks_queue;
    std::mutex queue_mutex;

    uint n_assemblies_maf = 0;
    std::thread parsing_thread(maf_parser, std::ref(maf_file), std::ref(blocks_queue), std::ref(queue_mutex), std::ref(parsing_ended), std::ref(n_assemblies_maf), false);

    std::thread processing_thread(processor, std::ref(blocks_queue), std::ref(queue_mutex), std::ref(ref_coverage), std::ref(parsing_ended), 100);

    parsing_thread.join();
    processing_thread.join();
    maf_file.close();

    std::unordered_map<std::string, std::unordered_map<uint, uint>> distribution;
    for (auto& scaffold: ref_coverage) {
        for (auto& cov: scaffold.second) {
            ++distribution[scaffold.first][cov];
        }
    }

    std::cout << "Contig\tCoverage\tCount\n";
    for (auto& scaffold: distribution) {
        for (auto& cov: scaffold.second) {
            std::cout << scaffold.first << "\t" << cov.first << "\t" << cov.second << "\n";
        }
    }
}

