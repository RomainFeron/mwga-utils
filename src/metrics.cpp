#include "metrics.h"


static const char USAGE[] =
R"(Generate wig files with base metrics from a MAF file.

    Usage:
      metrics <maf_file> [-p <prefix> -t <threads> -n <assemblies>]

    Options:
      <maf_file>       Path to a MAF file.
      -p <prefix>      Prefix for output wig files [default: no prefix]
      -n <assemblies>  Manually specify the number of assemblies in the alignment; if not, it is computed from the MAF [default: 0]
      -t <threads>     Number of threads to use [default: 1].
      -h --help        Show this screen.
)";


void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, std::mutex& results_mutex, bool& parsing_ended, Metrics& metrics, unsigned long batch_size) {

    uint pos = 0;
    size_t sep_pos = 0;
    std::string scaffold = "";
    std::vector<uint> tmp_alignability, tmp_identity;
    std::vector<MafBlock> batch;
    bool keep_going = true;

    while (keep_going) {
        batch = get_batch(blocks_queue, queue_mutex, batch_size);
        if (batch.size() > 0) {
            for (auto block: batch) {
                sep_pos = block.records[0].source.find(".");
                scaffold = block.records[0].source.substr(sep_pos + 1);
                tmp_alignability.resize(0);
                tmp_alignability.resize(block.records[0].length);
                tmp_identity.resize(0);
                tmp_identity.resize(block.records[0].length);
                pos = 0;
                for (uint i=0; i<block.records[0].sequence.size(); ++i) {  // Iterate over all positions in ref assembly for this block
                    if (block.records[0].sequence[i] != '-') {  // Exclude positions with gaps in ref assembly
                        for (auto record: block.records) {  // Iterate over all assemblies
                            if (record.sequence[i] != '-') ++tmp_alignability[pos];  // If no gap, position was aligned in this assembly
                            if (record.sequence[i] == block.records[0].sequence[i]) ++tmp_identity[pos];  // Check if ref and non-ref have the same nucleotide at this position
                        }
                        ++pos;
                    }
                }
                results_mutex.lock();
                for (uint i=0; i<block.records[0].length; ++i) {
                    if (metrics.find(scaffold) == metrics.end()) {
                        metrics[scaffold].resize(block.records[0].source_length);
                    }
                    metrics[scaffold][block.records[0].start + i].alignability += tmp_alignability[i];
                    metrics[scaffold][block.records[0].start + i].identity += tmp_identity[i];
                }
                results_mutex.unlock();
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
        std::cerr << "Error opening wig file <" << output_file_path << ">." << std::endl;
        exit(1);
    }
    // Format floating point in output stream to show 3 significant digits
    output_file << std::fixed << std::showpoint;
    output_file << std::setprecision(3);
    return output_file;
}



int main(int argc, char *argv[]) {

    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "Parser 0.1");
    std::string maf_file_path = args["<maf_file>"].asString();
    std::string prefix = args["-p"].asString();
    // Set empty prefix as default
    if (prefix == "no prefix") {
        prefix = "";
    }
    uint n_threads = static_cast<uint>(std::stoi(args["-t"].asString()));
    uint n_assemblies = static_cast<uint>(std::stoi(args["-n"].asString()));

    std::string alignability_wig_path = prefix + "alignability.wig";
    std::string identity_wig_path = prefix + "identity.wig";

    // Open input file
    std::ifstream maf_file = check_open(maf_file_path);
    bool parsing_ended = false;

    BlocksQueue blocks_queue;
    std::mutex queue_mutex, results_mutex;

    uint n_assemblies_maf = 0;
    std::thread parsing_thread(maf_parser, std::ref(maf_file), std::ref(blocks_queue), std::ref(queue_mutex), std::ref(parsing_ended), std::ref(n_assemblies_maf), false);

    Metrics metrics;

    std::vector<std::thread> processing_threads;
    for (uint i=0; i < n_threads; ++i) {
        processing_threads.push_back(std::thread(processor, std::ref(blocks_queue), std::ref(queue_mutex), std::ref(results_mutex), std::ref(parsing_ended), std::ref(metrics), 100));
    }

    parsing_thread.join();
    for (auto &t: processing_threads) t.join();

    if (n_assemblies == 0) n_assemblies = n_assemblies_maf;

    // Generate wig files
    std::cerr << "Generating wig files ..." << std::endl;
    std::ofstream alignability_wig = open_ofile(alignability_wig_path);
    std::ofstream identity_wig = open_ofile(identity_wig_path);
    std::string step_header = "";
    for (auto& scaffold: metrics) {
        // Write header for new scaffold
        step_header = "fixedStep chrom=" + scaffold.first + " start=1 step=1\n";
        alignability_wig << step_header;
        identity_wig << step_header;
        // Write values for each position in the scaffold
        for (auto& pos: scaffold.second) {
            alignability_wig << static_cast<float>(pos.alignability) / n_assemblies << "\n";
            identity_wig << static_cast<float>(pos.identity) / n_assemblies << "\n";
        }
    }

}
