#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "output.h"
#include "parser.h"
#include "utils.h"


typedef std::unordered_map<std::string, std::vector<BaseData>> MafMetrics;

void run_mafstats(int argc, char* argv[]);
void stats(Parameters& parameters);

typedef void (*command)(Parameters& parameters);
static std::unordered_map<std::string, command> commands {{"stats", &stats}};

//class MafStats {

//    public:

//        typedef void (MafStats::*command)();
//        std::unordered_map<std::string, command> commands {{"stats", &MafStats::stats}};

//        Parameters parameters;
//        std::unordered_map<std::string, std::vector<BaseData>> metrics;  // Store alignability results

//        MafStats();
//        MafStats(int argc, char* argv[]);
//        void stats();
//        void generate_output();
//};

