#include "stats.h"

static const char USAGE[] =
R"(Compute a series of statistics on a MAF file:
        - Number of BP aligned in each assembly

    Usage:
      stats <maf_file> [-p <prefix>]

    Options:
      <maf_file>       Path to a MAF file.
      -p <prefix>      Prefix for output stats files [default: stats]
      -h --help        Show this screen.
)";


void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, BpAlignedData& bp_aligned, bool& parsing_ended, ulong batch_size) {

    size_t sep_pos = 0;
    std::string scaffold = "", assembly = "";
    std::vector<MafBlock> batch;
    bool keep_going = true;

    while (keep_going) {
        batch = get_batch(blocks_queue, queue_mutex, batch_size);
        if (batch.size() > 0) {
            for (auto block: batch) {
                for (auto record: block.records) {
                    sep_pos = record.source.find(".");
                    assembly = record.source.substr(0, sep_pos);
                    if (bp_aligned.find(assembly) == bp_aligned.end()) bp_aligned[assembly] = 0;
                    bp_aligned[assembly] += record.length;
                }
            }
        }
        queue_mutex.lock();
        if (parsing_ended and blocks_queue.size() == 0) keep_going = false;
        queue_mutex.unlock();
    }
}


std::ofstream open_ofile(std::string& output_file_path) {
    std::ofstream output_file;
    output_file.open(output_file_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening stats file <" << output_file_path << ">." << std::endl;
        exit(1);
    }
    return output_file;
}



int main(int argc, char *argv[]) {

    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "Parser 0.1");
    std::string maf_file_path = args["<maf_file>"].asString();
    std::string prefix = args["-p"].asString();

    std::string bp_aligned_path = prefix + "_bp_aligned.tsv";

    BpAlignedData bp_aligned;

    // Process MAF file
    std::cerr << "Processing the MAF file" << std::endl;
    std::ifstream maf_file = check_open(maf_file_path);
    bool parsing_ended = false;

    BlocksQueue blocks_queue;
    std::mutex queue_mutex;

    uint n_assemblies_maf = 0;
    std::thread parsing_thread(maf_parser, std::ref(maf_file), std::ref(blocks_queue), std::ref(queue_mutex), std::ref(parsing_ended), std::ref(n_assemblies_maf), false);
    std::thread processing_thread(processor, std::ref(blocks_queue), std::ref(queue_mutex), std::ref(bp_aligned), std::ref(parsing_ended), 100);

    parsing_thread.join();
    processing_thread.join();
    maf_file.close();

    std::ofstream bp_aligned_file = open_ofile(bp_aligned_path);
    bp_aligned_file << "Assembly\tBp_aligned\n";
    for (auto& assembly: bp_aligned) {
        bp_aligned_file << assembly.first << "\t" << assembly.second << "\n";
    }
    bp_aligned_file.close();
}

