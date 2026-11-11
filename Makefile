CXX      ?= g++
CXXFLAGS ?= -O2 -pthread

BIN := bin/matmul

all: $(BIN)

$(BIN): src/matmul.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@

# Quick correctness checks on small sizes, including non-dividing thread counts
test: $(BIN)
	./$(BIN) 64 64 64 0 0 0 3
	./$(BIN) 64 64 64 0 0 16 3
	./$(BIN) 100 70 130 7 1 16 3
	./$(BIN) 100 70 130 3 4 0 3
	./$(BIN) 1000 128 200 6 6 32 3
	./$(BIN) 50 50 50 36 1 8 3

version:
	$(CXX) --version | head -1

clean:
	rm -rf build

.PHONY: all test version clean