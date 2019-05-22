#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "coverage.h"
#include "output.h"
#include "parser.h"
#include "stats.h"
#include "utils.h"


void run_mafstats(int argc, char* argv[]);

typedef void (*command)(Parameters& parameters);
static std::unordered_map<std::string, command> commands {{"stats", &stats},
                                                          {"cov", &coverage}};

