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
TARGETS = metrics

# Rules
all: init $(TARGETS)

test: $(SRC)/test.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/test $^ $(LDFLAGS)

metrics: $(SRC)/metrics.cpp $(SRC)/maf_parser.cpp $(INCLUDE)/docopt/docopt.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/metrics $^ $(LDFLAGS)

# $(TARGET): $(OBJS) $(LIBOBJS)
# 	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/$(TARGET) $^ $(LDFLAGS)

# $(BUILD)/%.o: $(SRC)/%.cpp
# 	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

# $(LIBBUILD)/%.o: $(INCLUDE)/*/%.cpp
# 	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

clean:
	@rm -rf $(BUILD)/*.o
	@rm -rf $(BIN)/$(TARGETS)
	@rm -rf $(INCLUDE)/*/*.o

init:
	@mkdir -p $(BUILD) $(BUILD)
	@mkdir -p $(BIN) $(BIN)

rebuild: clean $(TARGETS)
