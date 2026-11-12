CXX ?= g++

# Directory containing assignment2.h (starter kit); override with: make INC=path
INC ?= src/partA

PTHREADS := build/matmul_pthreads
OMP      := build/matmul_omp

all: pthreads omp

pthreads: $(PTHREADS)
omp: $(OMP)

$(PTHREADS): src/partA/matmul.cpp
	@mkdir -p build
	$(CXX) -O2 -pthread -I$(INC) $< -o $@

$(OMP): src/partB/matmul_omp.cpp
	@mkdir -p build
	$(CXX) -O2 -fopenmp -I$(INC) $< -o $@

# Quick correctness checks on small sizes, including non-dividing thread counts
test: all
	./$(PTHREADS) 64 64 64 0 0 0 3
	./$(PTHREADS) 100 70 130 7 1 16 3
	./$(PTHREADS) 1000 128 200 6 6 32 3
	OMP_PROC_BIND=true OMP_PLACES=cores ./$(OMP) 64 64 64 0 0 0 3
	OMP_PROC_BIND=true OMP_PLACES=cores ./$(OMP) 100 70 130 7 1 16 3
	OMP_PROC_BIND=true OMP_PLACES=cores ./$(OMP) 1000 128 200 6 6 32 3
	OMP_PROC_BIND=true OMP_PLACES=cores ./$(OMP) 50 50 50 36 1 8 3

version:
	$(CXX) --version | head -1

clean:
	rm -rf build

.PHONY: all pthreads omp test version clean