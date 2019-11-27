#include "utils.h"

// Output current date and time in format specified in utils.h
char* print_time (char *buff) {

    time_t t = time (nullptr);
    strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
    return buff;
}



// Splits a std::string into a std::vector of std::strings according to a specified delimiter (default: \t)
std::vector<std::string> split(std::string str, const std::string delimiter){

    std::vector<std::string> output;
    size_t pos;

    while ((pos = str.find(delimiter)) != std::string::npos){

        output.push_back(str.substr(0, pos));
        str.erase(0, pos + delimiter.length());
    }

    output.push_back(str.substr(0, pos));

    return output;
}


void rev_comp(const std::string& sequence, std::string& revcomp_sequence) {

    revcomp_sequence = "";
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        revcomp_sequence += revcomp_table[int(*it)];
    }
}


std::ifstream check_open(std::string input_file_path) {

    std::ifstream input_file;
    input_file.open(input_file_path);

    if (not input_file.is_open()) {
        std::cerr << "Error opening input file <" << input_file_path << ">." << std::endl;
        exit(1);
    }

    return input_file;
}


// Faster string to int conversion
inline int fast_stoi(const char* str) {

    int val = 0;
    while( *str ) {
        val = val*10 + (*str++ - '0');
    }
    return val;
}
