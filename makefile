# Compiler options
CC = g++
OPTCFLAGS = -O2
CFLAGS = -Wall -std=gnu++11 $(OPTCFLAGS)
LDFLAGS = -pthread -static-libstdc++ -lz

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
BUILD = $(BASEDIR)/build
LIBBUILD = $(BASEDIR)/build
INCLUDE = $(BASEDIR)/include
SRC = $(BASEDIR)/src
CPP = $(wildcard $(SRC)/*.cpp) $(wildcard $(SRC)/*/*.cpp)
LIBCPP = $(wildcard $(INCLUDE)/*/*.cpp)

# Targets
TARGETS = metrics missing_regions

# Rules
all: init $(TARGETS)

metrics: $(SRC)/metrics.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/metrics $^ $(LDFLAGS)

missing_regions: $(SRC)/missing_regions.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/missing_regions $^ $(LDFLAGS)

clean:
	@rm -rf $(BUILD)/*.o
	@rm -rf $(BIN)/$(TARGETS)
	@rm -rf $(INCLUDE)/*/*.o

init:
	@mkdir -p $(BUILD) $(BUILD)
	@mkdir -p $(BIN) $(BIN)

rebuild: clean $(TARGETS)
