#include "missing_regions.h"


static const char USAGE[] =
R"(Add regions from the reference genome that are missing from a MAF file.

    Usage:
      missing_regions <maf_file> <reference>

    Options:
      <maf_file>       Path to a MAF file.
      <reference>      Path to a FASTA file for the reference assembly.
      -h --help        Show this screen.
)";


void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, MafCoverage& maf_coverage, bool& parsing_ended, ulong batch_size) {

    size_t sep_pos = 0;
    std::string scaffold = "";
    std::vector<MafBlock> batch;
    bool keep_going = true;
    MafCoverage tmp_maf_coverage;

    while (keep_going) {
        batch = get_batch(blocks_queue, queue_mutex, batch_size);
        if (batch.size() > 0) {
            for (auto block: batch) {
                sep_pos = block.records[0].source.find(".");
                scaffold = block.records[0].source.substr(sep_pos + 1);
                if (maf_coverage.find(scaffold) == maf_coverage.end()) {
                    maf_coverage[scaffold].resize(block.records[0].source_length);
                }
                for (uint i=0; i<block.records[0].length; ++i) {
                    maf_coverage[scaffold][block.records[0].start + i] = 1;
                }
                std::cout << "a score=" << block.score << "\n";
                for (auto record: block.records) std::cout << "s " << record.source << " " << record.start << " " << record.length << " " << record.strand << " "
                                                            << record.source_length << " " << record.sequence << "\n";
                std::cout << "\n";
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
    std::string reference_file_path = args["<reference>"].asString();

    MafCoverage maf_coverage;

    // Process genome file
    std::cerr << "Processing the reference assembly" << std::endl;
    std::ifstream reference_file = check_open(reference_file_path);
    Assembly reference;
    std::string line = "", header = "", sequence = "";
    std::vector<std::string> tmp;
    while (std::getline(reference_file, line)) {
        if (line[0] == '>') {
            if (header != "") {
                reference[header] = sequence;
                maf_coverage[header] = std::vector<bool>(sequence.size(), 0);
            }
            header = line.substr(1);
            tmp = split(header, ":");
            header = tmp[1];
            sequence = "";
        } else {
            sequence += line;
        }
    }
    reference[header] = sequence;
    reference_file.close();

    // Process MAF file
    std::cerr << "Processing the MAF file" << std::endl;
    std::ifstream maf_file = check_open(maf_file_path);
    bool parsing_ended = false;

    BlocksQueue blocks_queue;
    std::mutex queue_mutex;

    uint n_assemblies_maf = 0;
    std::thread parsing_thread(maf_parser, std::ref(maf_file), std::ref(blocks_queue), std::ref(queue_mutex), std::ref(parsing_ended), std::ref(n_assemblies_maf));
    std::thread processing_thread(processor, std::ref(blocks_queue), std::ref(queue_mutex), std::ref(maf_coverage), std::ref(parsing_ended), 100);

    parsing_thread.join();
    processing_thread.join();

    uint start = 0, length = 0, scaffold_size = 0;
    std::string scaffold_name = "";
    bool new_sequence = true;

    std::cerr << "Exporting regions absent from the MAF file" << std::endl;

    for (auto scaffold: maf_coverage) {
        sequence = "";
        start = 0;
        length = 0;
        new_sequence = true;
        scaffold_name = scaffold.first;
        scaffold_size = static_cast<uint>(scaffold.second.size());
        for (uint i=0; i<scaffold_size; ++i) {
            if (scaffold.second[i] or i == scaffold.second.size() - 1) {
                if (not new_sequence) {
                    sequence = reference[scaffold_name].substr(start, length);
                    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                    std::cout << "a score=NA\n";
                    std::cout << "s " << scaffold_name << " " << start << " " << length << " + " << scaffold_size << " " << sequence << "\n\n";
                }
                sequence = "";
                start = 0;
                length = 0;
                new_sequence = true;
            } else {
                if (new_sequence) {
                    start = i;
                    new_sequence = false;
                }
                ++length;
            }
        }
    }

    maf_file.close();
}
