# Multilevel DAG Partitioning

A C++ implementation of a multilevel algorithm for partitioning Directed Acyclic Graphs (DAGs). This work is based on the paper ["Multilevel Algorithms for Acyclic Partitioning of Directed Acyclic Graphs"](https://epubs.siam.org/doi/abs/10.1137/18M1176865) by Herrmann et al.

The tool focuses on producing high-quality partitions by minimizing the edge cut and balancing partition loads, all while preserving the critical acyclic dependencies between partitions.

## How It Works

The algorithm works in three main phases to efficiently partition large DAGs:

1. **Coarsening:** The original graph is progressively simplified into smaller, approximate graphs by clustering nodes, preserving the overall structure and acyclicity.
2. **Initial Partitioning:** The smallest graph is partitioned into two parts.
3. **Refinement:** The partition is then projected back to the original, larger graph, with refinement heuristics applied at each step to improve the partition quality.

This process uses recursive bisection to achieve partitions into `k` parts.

## Features

- **High-performance containers**: Uses Robin Hood hashing for improved memory efficiency and performance
- **OpenMP support**: Optional parallel processing capabilities
- **Memory-aware scheduling**: After partitioning, optionally schedule the execution order of partitions to minimize peak memory usage using a CP-SAT solver

## Requirements

- A C++17 compatible compiler.
- CMake 3.28 or later.
- [OR-Tools](https://developers.google.com/optimization) (system installation required, used for the CP-SAT scheduler).
- (Optional) Valgrind for memory profiling.

## Building the Project

This project uses CMake to automatically download and build most of its dependencies (GKlib, METIS, Scotch, and Robin Hood Hashing). OR-Tools must be installed separately as a system dependency before building.

### Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```

### Installation

To install the library and executables to a custom path:

```bash
cmake -DCMAKE_INSTALL_PREFIX=../install ..
make
make install
```

This will install the library, headers, and executables to the `install` directory.

### Build Options

- `DAG_PARTITIONING_OPENMP=ON` - Enable OpenMP support for parallel processing
- `BUILD_SHARED_LIBS=ON` - Build shared libraries instead of static
- `CMAKE_BUILD_TYPE=Release` - Build optimized release version (default)

Example with options (for more options see `CMakeLists.txt`):

```bash
cmake -DDAG_PARTITIONING_OPENMP=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=../install ..
make
make install
```

## Usage

### Using the Test Executable

After building, you can use the `dag-test` executable to partition DAGs:

```bash
# From the install directory

# Required arguments:
#   <# of partitions>       Number of partitions to produce
#   <dot file path>         Path to the input DOT file
#   <clustering method>     FORB, CYC, or HYB
#   <bisection method>      GGG, UNDIRSCOTCH, UNDIRMETIS, or UNDIRBOTH
#   <refinement method>     BOUNDARYFM, BOUNDARYKL, or MIXED
#   <enable parallel>       1 to enable parallel processing, 0 for sequential
# Optional arguments:
#   [min size for parallel] Minimum subgraph size to parallelize (default: 100)
#   [max parallel depth]    Maximum recursion depth for parallelization (default: 10)
#   [enable scheduling]     1 to run the memory-aware scheduler after partitioning, 0 to skip (default: 0)

# Example: partition into 4 parts with default settings, parallel enabled
./bin/dag-test 4 ../test/example.dot HYB GGG BOUNDARYFM 1

# Example: also run the scheduler
./bin/dag-test 4 ../test/example.dot HYB GGG BOUNDARYFM 1 100 10 1

# For detailed usage
./bin/dag-test
```

### Generating Random DAGs

Use the `rand-dag` tool to generate test graphs:

```bash
# From the install directory

# Generate a random DAG with 100 nodes, 150% edge density and 0 debug level
./bin/rand-dag 1000 150 0

# This creates a DOT format graph that can be used with the dag-test executable
./bin/rand-dag 1000 150 0 random_dag.dot
./bin/dag-test 4 random_dag.dot HYB GGG BOUNDARYFM 1

# For detailed usage
./bin/rand-dag
```

### Library Usage

For programmatic usage, check the dag-test executable source in [test/test.cpp](test/test.cpp).

## Input Format

The input DAG should be provided in the dot format as follows:

- First, the lines for the vertices: *vertexname*[weight=*vertexweight*];
- Then, lines for the edges: *from*->*to*[weight=*edgeweight*];

Example:

```
digraph cfg {
0[weight=1];
1[weight=1];
2[weight=1];
3[weight=1];
4[weight=1];
5[weight=1];
6[weight=1];
7[weight=1];
8[weight=1];
9[weight=1];
10[weight=1];
0->1[weight=1];
0->6[weight=1];
0->7[weight=1];
1->2[weight=1];
1->5[weight=1];
1->7[weight=1];
1->8[weight=1];
1->10[weight=1];
3->9[weight=1];
4->7[weight=1];
5->8[weight=1];
6->3[weight=1];
6->9[weight=1];
8->2[weight=1];
8->9[weight=1];
9->4[weight=1];
}
```

## Scheduling

After partitioning, the tool can optionally compute an execution order for the partitions that minimizes peak memory usage. This is useful when the partitions will be executed sequentially and memory is a constraint.

The scheduler works as follows:

1. **Coarse graph construction**: Builds a graph where each node represents a partition. Node weights capture the peak memory used within each partition (computed via best-fit tensor packing), and edge weights represent the volume of inter-partition data transfers.
2. **CP-SAT formulation**: The partition execution order is found by solving a constraint programming problem using Google OR-Tools' CP-SAT solver. The formulation is an adaptation to CP-SAT of an ILP formulation from [Zhong et.al.](https://arxiv.org/abs/2308.13898).
3. **Warm start**: An RPO (Reverse Post-Order) heuristic provides an initial feasible schedule that seeds the CP-SAT solver, accelerating convergence.

The scheduler is invoked by passing `1` as the last argument to `dag-test` (see Usage above). It reports the peak memory of the found schedule and the time taken by the solver.

## Output Format

The program outputs (comma-separated):

- Clustering, bisection, and refinement method names
- Edge cut value
- Partition load imbalance (%)
- Partitioning time (microseconds)
- Scheduling time in microseconds (0 if scheduling was not run)

## Performance and Comparison

This implementation was benchmarked against the original dagP implementation. The key findings are:

- Stability: This implementation proved to be more stable, completing all tests without the random crashes encountered with the original tool.

- Performance: With multi-threading enabled, this implementation is significantly faster for many workloads.

- Quality: Partition quality (edge cut and load balance) is competitive with the original.

**A Note on Fair Comparison:** The original dagP tool is single-threaded. Our benchmarks run this implementation with multi-threading enabled by default, giving it a performance advantage. For a pure single-threaded comparison, you can disable multi-threading when running dag-test.

### Reproducing the Benchmarks

Scripts to reproduce our performance and memory profiling results are in the test/ directory. You will need to build the original dagP library for a full comparison.

- See the [`run-compare.py`](test/run-compare.py) script for performance comparison.
- See the [`run-memprof.py`](test/run-memprof.py) script for memory profiling with Valgrind.

For more information about the benchmarks, see [`TEST.md`](TEST.md).

## License

CC BY-NC 4.0

## References

Partitioning based on the [paper](https://epubs.siam.org/doi/abs/10.1137/18M1176865):
"Multilevel Algorithms for Acyclic Partitioning of Directed Acyclic Graphs" by Herrmann et al.

Scheduling formulation adapted to CP-SAT from [Zhong et.al.](https://arxiv.org/abs/2308.13898).
