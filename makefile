# Compiler options
CC = g++
OPTCFLAGS = -O2
CFLAGS = -Wall -std=c++11 $(OPTCFLAGS)
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

# Target
TARGET = maf_stats

# Variables
OBJS = $(addprefix $(BUILD)/, $(notdir $(CPP:.cpp=.o)))
LIBOBJS = $(LIBCPP:.cpp=.o)

# Rules

all: init $(TARGET)

print-%: ; @echo $* = $($*)

$(TARGET): $(OBJS) $(LIBOBJS)
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/$(TARGET) $^ $(LDFLAGS)

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

$(LIBBUILD)/%.o: $(INCLUDE)/*/%.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

clean:
	@rm -rf $(BUILD)/*.o
	@rm -rf $(BIN)/$(TARGET)
	@rm -rf $(INCLUDE)/*/*.o

init:
	@mkdir -p $(BUILD) $(BUILD)
	@mkdir -p $(BIN) $(BIN)

rebuild: clean $(TARGET)
