#include <thread>
#include "metrics.h"
#include "maf_parser.h"


static const char USAGE[] =
R"(Generate wig files with base metrics from a MAF file.

    Usage:
      metrics <maf_file> [-p <prefix> -t <threads>]

    Options:
      <maf_file>    Path to a MAF file.
      -p <prefix>   Prefix for output wig files [default: metrics]
      -t <threads>  Number of threads to use [default: 1].
      -h --help     Show this screen.
)";


void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, std::mutex& results_mutex, bool& parsing_ended, Metric& alignability, ulong batch_size = 100) {

    uint pos = 0;
    size_t sep_pos = 0;
    std::string scaffold = "";
    std::vector<uint> tmp_metric;
    std::vector<MafBlock> batch;
    bool keep_going = true;

    while (keep_going) {
        batch = get_batch(blocks_queue, queue_mutex, batch_size);
        if (batch.size() > 0) {
            for (auto block: batch) {
                sep_pos = block.records[0].source.find(".");
                scaffold = block.records[0].source.substr(sep_pos + 1);
                tmp_metric.resize(0);
                tmp_metric.resize(block.records[0].length);
                pos = 0;
                for (uint i=0; i<block.records[0].sequence.size(); ++i) {
                    if (block.records[0].sequence[i] != '-') {
                        for (auto record: block.records) {
                            if (record.sequence[i] != '-') ++tmp_metric[pos];
                        }
                        ++pos;
                    }
                }
                results_mutex.lock();
                for (uint i=0; i<block.records[0].length; ++i) {
                    if (alignability.find(scaffold) == alignability.end()) {
                        alignability[scaffold].resize(block.records[0].source_length);
                    }
                    alignability[scaffold][block.records[0].start + i] += tmp_metric[i];
                }
                results_mutex.unlock();
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
    std::string prefix = args["-p"].asString();
    uint n_threads = static_cast<uint>(std::stoi(args["-t"].asString()));

    std::string alignability_wig_path = prefix + "_alignability.wig";

    // Open input file
    std::ifstream maf_file = check_open(maf_file_path);
    bool parsing_ended = false;

    BlocksQueue blocks_queue;
    std::mutex queue_mutex, results_mutex;

    std::thread parsing_thread(maf_parser, std::ref(maf_file), std::ref(blocks_queue), std::ref(queue_mutex), std::ref(parsing_ended));

    Metric alignability;

    std::vector<std::thread> processing_threads;
    for (uint i=0; i < n_threads; ++i) {
        processing_threads.push_back(std::thread(processor, std::ref(blocks_queue), std::ref(queue_mutex), std::ref(results_mutex), std::ref(parsing_ended), std::ref(alignability), 100));
    }

    parsing_thread.join();
    for (auto &t: processing_threads) t.join();

    std::cerr << "Generating alignability wig file ..." << std::endl;

    // Open wig file
    std::ofstream output_file;
    output_file.open(alignability_wig_path);
    if (not output_file.is_open()) {
        std::cerr << "Error opening wig file <" << alignability_wig_path << ">." << std::endl;
        exit(1);
    }

    // Output the data
    for (auto& scaffold: alignability) {
        output_file << "fixedStep chrom=" << scaffold.first << " start=1 step=1\n";
        for (auto& pos: scaffold.second) {
            output_file << pos << "\n";
        }
    }

}
