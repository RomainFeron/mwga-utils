#include "maf_parser.h"

void maf_parser(std::ifstream& maf_file, BlocksQueue& blocks_queue, std::mutex& queue_mutex, bool& parsing_ended, uint& n_assemblies, bool print_comments) {
    /* Parse a MAF file and store alignment blocks in a queue shared by processing threads.
     * Each sequence in a block is stored in a MafRecord struct, and each block is stored in a MafBlock struct.
     * Blocks are first stored in a temporary queue which is merged with the shared queue every 1000 blocks.
     */

    char buffer[65536];  // Input buffer where file is read and stored
    long k = 0;  // Real size of the block that was read (for iterating over the buffer)
    char previous_char = ' ';  // Store the last character processed from the buffer
    std::string temp = "", comment = "";  // Store the current field string
    uint field_n = 0, blocks_n = 0;  // Fields and blocks counters
    bool new_line = true, new_block = false;  // Flags for new line and new block
    MafRecord maf_record;  // Structure storing record data (i.e. single sequence within block)
    MafBlock maf_block;  // Structure storing block data
    std::vector<MafBlock> tmp_queue(1000);  // Temporary block queue to avoid locking the shared blocks queue too often

    n_assemblies = 0;  // Number of assemblies in the MAF file
    bool in_comment = false;

    do {

        // Read a chunk of the maf file into a buffer and get the real length of this chunk
        maf_file.read(buffer, sizeof(buffer));
        k = maf_file.gcount();

        for (uint i=0; i<k; ++i) {  // Iterate over data in buffer

            // Process a single character
            switch(buffer[i]) {

                case ' ':  // Fields are delimited by varying number of spaces
                    if (previous_char != ' ') {  // If previous character was not a space, we reached a delimiter and we save the field to the current record
                        switch (field_n) {
                            case 1:
                                maf_record.source = temp;  // Second field is the source (assembly + contig)
                                break;
                            case 2:
                                maf_record.start = static_cast<uint>(fast_stoi(temp.c_str()));  // Third field is the start of this block in this assembly
                                break;
                            case 3:
                                maf_record.length = static_cast<uint>(fast_stoi(temp.c_str()));  // Third field is the length of this block in this assembly
                                break;
                            case 4:
                                maf_record.strand = temp[0];  // Fourth field is the strand of this block in this assembly
                                break;
                            case 5:
                                maf_record.source_length = static_cast<uint>(fast_stoi(temp.c_str()));  // Fifth field is the length of the contig containing this block in this assembly
                                break;
                            default:
                                break;
                        }
                        ++field_n;
                        temp = "";  // Reset string storing field data as we finished processing a field
                    }
                    if (in_comment) comment += buffer[i];
                    break;

                case '\n':

                    if (new_block) {
                        maf_block.score = temp.substr(6); // Save block score
                        new_block = false;
                    }

                    if (field_n == 6) { // This line was a MAF record (there are only 7 fields, 0-based indexing)
                        maf_record.sequence = temp;  // Add sequence to record
                        maf_block.records.push_back(maf_record); // Add record to current block
                        ++maf_block.n_records;
                    }

                    if (new_line) { // This line is empty (previous character was already '\n')
                        tmp_queue[blocks_n % 1000] = maf_block;  // Empty line means end of a block, we add it to the queue
                        if (maf_block.n_records > n_assemblies) n_assemblies = maf_block.n_records;  // Update max number of assemblies in a block
                        maf_block.records.resize(0);  // Reset record data
                        maf_block.n_records = 0;
                        ++blocks_n;
                        if (blocks_n % 1000 == 0) {  // Merge temporary queue with shared queue after 1000 blocks
                            if (blocks_n % 10000 == 0) std::cerr << "Processed " << blocks_n / 1000 << " K. blocks" << std::endl;
                            queue_mutex.lock();
                            for (auto block: tmp_queue) blocks_queue.push(block);
                            queue_mutex.unlock();
                            tmp_queue.resize(0);  // Reset temporary queue
                            tmp_queue.resize(1000);
                        }

                    }

                    new_line = true;

                    if (print_comments and in_comment and comment != "##eof maf") std::cout << comment << std::endl;  // Assumption: this is safe because all comments should be within the first 1000 lines of the MAF

                    // Reset variables
                    field_n = 0;
                    temp = "";
                    in_comment = false;
                    comment = "";
                    break;

                case '#':
                    in_comment = true;
                    comment += buffer[i];
                    break;

                default:

                    if (new_line) {
                        switch (buffer[i]) {
                            case 'a':  // Line starting with 'a' marks the beginning of a block
                                new_block = true;
                                break;
                            default:
                                break;
                        }
                        new_line = false;
                    } else {
                        temp += buffer[i];
                    }
                    if (in_comment) comment += buffer[i];
                    break;
            }

            previous_char = buffer[i];
        }

    } while (maf_file);

    // Add the final blocks to the shared queue
    queue_mutex.lock();
    for (auto block: tmp_queue) {
        if (block.n_records > 0) blocks_queue.push(block);
    }
    queue_mutex.unlock();

    parsing_ended = true;  // Shared flag indicating that the parsing is finished

}



std::vector<MafBlock> get_batch(BlocksQueue& blocks_queue, std::mutex& queue_mutex, unsigned long batch_size) {
    /* Get a batch of <batch_size> blocks from the shared queue.
     * The batch is stored as a vector of blocks. Extracted blocks are removed from the shared queue.
     */

    queue_mutex.lock();

    unsigned long batch_size_real = std::min(batch_size, blocks_queue.size());
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
