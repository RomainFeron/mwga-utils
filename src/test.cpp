#include "maf_parser.h"
#include <thread>

static const char USAGE[] =
R"(Test.

    Usage:
      test <maf_file> [-t <threads>]

    Options:
      <maf_file>    Path to a MAF file.
      -t <threads>  Number of threads to use [default: 1].
      -h --help     Show this screen.
)";


void processor(BlocksQueue& blocks_queue, std::mutex& queue_mutex, std::mutex& score_mutex, bool& parsing_ended, float& total_score, ulong batch_size = 100) {

    float tmp_score = 0;
    std::vector<MafBlock> batch;
    bool keep_going = true;
    while (keep_going) {
        tmp_score = 0;
        batch = get_batch(blocks_queue, queue_mutex, batch_size);
        if (batch.size() > 0) {
            for (auto block: batch) {
                tmp_score += block.score;
            }
        }
        score_mutex.lock();
        total_score += tmp_score;
        score_mutex.unlock();
        queue_mutex.lock();
        if (parsing_ended and blocks_queue.size() == 0) keep_going = false;
        queue_mutex.unlock();
    }
}


int main(int argc, char *argv[]) {

    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "Parser 0.1");
    std::string maf_file_path = args["<maf_file>"].asString();
    uint n_threads = static_cast<uint>(std::stoi(args["-t"].asString()));

    // Open input file
    std::ifstream maf_file = check_open(maf_file_path);
    float total_score = 0;
    bool parsing_ended = false;

    BlocksQueue blocks_queue;
    std::mutex queue_mutex, score_mutex;

    std::thread parsing_thread(maf_parser, std::ref(maf_file), std::ref(blocks_queue), std::ref(queue_mutex), std::ref(parsing_ended));

    std::vector<std::thread> processing_threads;
    for (uint i=0; i < n_threads; ++i) {
        processing_threads.push_back(std::thread(processor, std::ref(blocks_queue), std::ref(queue_mutex), std::ref(score_mutex), std::ref(parsing_ended), std::ref(total_score), 100));
    }

    parsing_thread.join();
    for (auto &t: processing_threads) t.join();

    std::cout << total_score << std::endl;

}
