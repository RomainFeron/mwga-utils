# Compiler options
CC = g++
OPTCFLAGS = -O2
CFLAGS += -Wall -std=c++11 $(OPTCFLAGS)
LDFLAGS += -pthread -lstdc++ -lm

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
BUILD = $(BASEDIR)/build
INCLUDE = $(BASEDIR)/include
SRC = $(BASEDIR)/src

# Targets
TARGETS = metrics missing_regions single_cov stats
# TARGETS = metrics

# Rules
all: init $(TARGETS)

$(INCLUDE)/docopt/docopt.o: $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $<

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $<

metrics: $(BUILD)/metrics.o $(BUILD)/maf_parser.o $(INCLUDE)/docopt/docopt.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/metrics $^ $(LDFLAGS)

missing_regions: $(BUILD)/missing_regions.o $(BUILD)/maf_parser.o $(INCLUDE)/docopt/docopt.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/missing_regions $^ $(LDFLAGS)

single_cov: $(BUILD)/single_cov.o $(BUILD)/maf_parser.o $(INCLUDE)/docopt/docopt.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/single_cov $^ $(LDFLAGS)

stats: $(BUILD)/stats.o $(BUILD)/maf_parser.o $(INCLUDE)/docopt/docopt.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/stats $^ $(LDFLAGS)

clean:
	@rm -rf $(BUILD)/*.o
	@rm -rf $(BIN)/$(TARGETS)
	@rm -rf $(INCLUDE)/*/*.o

init:
	@mkdir -p $(BUILD) $(BUILD)
	@mkdir -p $(BIN) $(BIN)

rebuild: clean $(TARGETS)
