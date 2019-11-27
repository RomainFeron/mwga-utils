#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <time.h>
#include <unordered_map>
#include <vector>
#define DTTMFMT "%Y-%m-%d %H:%M:%S"
#define DTTMSZ 21

const char revcomp_table[] = {
      0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
     16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
     32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
     48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
     64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

// Output current date and time in format specified in utils.h
inline char* print_time (char *buff) {

    time_t t = time (nullptr);
    strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
    return buff;
}


// Splits a std::string into a std::vector of std::strings according to a specified delimiter (default: \t)
inline std::vector<std::string> split(std::string str, const std::string delimiter){

    std::vector<std::string> output;
    size_t pos;

    while ((pos = str.find(delimiter)) != std::string::npos){

        output.push_back(str.substr(0, pos));
        str.erase(0, pos + delimiter.length());
    }

    output.push_back(str.substr(0, pos));

    return output;
}


inline void rev_comp(const std::string& sequence, std::string& revcomp_sequence) {

    revcomp_sequence = "";
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        revcomp_sequence += revcomp_table[int(*it)];
    }
}


inline std::ifstream check_open(std::string input_file_path) {

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


