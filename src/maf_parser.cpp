#include "maf_parser.h"

void maf_parser(std::ifstream& maf_file, BlocksQueue& blocks_queue, std::mutex& queue_mutex, bool& parsing_ended) {

    char buffer[65536];
    long k = 0;
    char previous_char = ' ';
    std::string temp = "";
    uint field_n = 0, blocks_n = 0;
    bool new_line = true, new_block = false;
    MafRecord maf_record;
    MafBlock maf_block;
    std::vector<MafBlock> tmp_queue(1000);

    do {

        maf_file.read(buffer, sizeof(buffer));
        k = maf_file.gcount();

        for (uint i=0; i<k; ++i) {

            // Read the buffer character by character
            switch(buffer[i]) {

                case ' ':  // Fields are delimited by varying number of spaces
                    if (previous_char != ' ') {  // If previous character was not a space, we save the field to the current record
                        switch (field_n) {
                            case 1:
                                maf_record.source = temp;
                                break;
                            case 2:
                                maf_record.start = static_cast<uint>(fast_stoi(temp.c_str()));
                                break;
                            case 3:
                                maf_record.length = static_cast<uint>(fast_stoi(temp.c_str()));
                                break;
                            case 4:
                                maf_record.strand = temp[0];
                                break;
                            case 5:
                                maf_record.source_length = static_cast<uint>(fast_stoi(temp.c_str()));
                                break;
                            default:
                                break;
                        }
                        ++field_n;
                        temp = "";
                    }
                    break;

                case '\n':

                    if (new_block) {
                        maf_block.score = std::stof(temp.substr(6)); // Save block score
                        new_block = false;
                    }

                    if (field_n == 6) { // This line was a MAF record
                        maf_record.sequence = temp;  // Add sequence to record
                        maf_block.records.push_back(maf_record); // Add record to current block
                        ++maf_block.n_records;
                    }

                    if (new_line) { // This line is empty (previous character was already '\n')
                        tmp_queue[blocks_n % 1000] = maf_block;  // Empty line means end of a block, we add it to the queue
                        maf_block.records.resize(0);
                        maf_block.n_records = 0;
                        ++blocks_n;
                        if (blocks_n % 1000 == 0) {
                            std::cerr << "Processed " << blocks_n / 1000 << " K. blocks" << std::endl;
                            queue_mutex.lock();
                            for (auto block: tmp_queue) blocks_queue.push(block);
                            queue_mutex.unlock();
                            tmp_queue.resize(0);
                            tmp_queue.resize(1000);
                        }
                    }

                    new_line = true;

                    // Reset variables
                    field_n = 0;
                    temp = "";
                    break;

                default:

                    if (new_line) {
                        switch (buffer[i]) {
                            case 'a':
                                new_block = true;
                                break;
                            default:
                                break;
                        }
                        new_line = false;
                    } else {
                        temp += buffer[i];
                    }

                    break;
            }

            previous_char = buffer[i];
        }

    } while (maf_file);

    queue_mutex.lock();
    for (auto block: tmp_queue) {
        if (block.n_records > 0) blocks_queue.push(block);
    }
    queue_mutex.unlock();

    parsing_ended = true;
}



std::vector<MafBlock> get_batch(BlocksQueue& blocks_queue, std::mutex& queue_mutex, ulong batch_size) {

    queue_mutex.lock();

    ulong batch_size_real = std::min(batch_size, blocks_queue.size());
    std::vector<MafBlock> batch(batch_size_real);

    if (batch_size_real > 0) {
        for (uint i=0; i<batch_size_real; ++i) {
            batch[i] = blocks_queue.front();
            blocks_queue.pop();
        }
    }

    queue_mutex.unlock();

    return batch;

}
