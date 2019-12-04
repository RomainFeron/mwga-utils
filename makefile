# Compiler options
CC = g++
OPTCFLAGS = -O2
CFLAGS = -Wall -std=gnu++11 $(OPTCFLAGS)
LDFLAGS = -pthread -lstdc++ -lm

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
BUILD = $(BASEDIR)/build
INCLUDE = $(BASEDIR)/include
SRC = $(BASEDIR)/src

# Targets
TARGETS = metrics missing_regions single_cov stats

# Rules
all: init $(TARGETS)

metrics: $(SRC)/metrics.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/metrics $^ $(LDFLAGS)

missing_regions: $(SRC)/missing_regions.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/missing_regions $^ $(LDFLAGS)

single_cov: $(SRC)/single_cov.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/single_coverage $^ $(LDFLAGS)

stats: $(SRC)/stats.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/stats $^ $(LDFLAGS)

clean:
	@rm -rf $(BUILD)/*.o
	@rm -rf $(BIN)/$(TARGETS)
	@rm -rf $(INCLUDE)/*/*.o

init:
	@mkdir -p $(BUILD) $(BUILD)
	@mkdir -p $(BIN) $(BIN)

rebuild: clean $(TARGETS)
