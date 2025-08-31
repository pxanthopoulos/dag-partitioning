# Multilevel DAG Partitioning Implementation

This repository contains a C++ implementation of a multilevel algorithm for partitioning Directed Acyclic Graphs (DAGs)
based on the research [paper](https://epubs.siam.org/doi/abs/10.1137/18M1176865) by Herrmann et al. The implementation
focuses on minimizing edge cuts and balancing the loads of the partitions, while maintaining acyclic dependencies
between partitions.

## Algorithm Overview

The algorithm implements a multilevel approach for partitioning DAGs with three main phases:

1. **Coarsening/Clustering Phase**:
    - Reduces the size of the input graph while preserving its essential structure
    - Uses specialized coarsening heuristics to maintain acyclicity
    - Creates a hierarchy of increasingly smaller graphs

2. **Initial Bisectioning Phase**:
    - Partitions the coarsest graph into 2 parts
    - Ensures balanced partition sizes
    - Maintains acyclic dependencies between partitions

3. **Refinement Phase**:
    - Projects the partition back through the hierarchy
    - Refines the partition at each level
    - Uses novel heuristics to improve the partition quality while preserving acyclicity

A recursive bisectioning scheme is employed to partition the
graph into the required number of parts.

## Features

- **High-performance containers**: Uses Robin Hood hashing for improved memory efficiency and performance
- **OpenMP support**: Optional parallel processing capabilities

## Requirements

- C++17 or later
- CMake 3.28 or later
- Git
- Valgrind for memory footprint profiling

## Building

This project automatically fetches and builds all required dependencies using CMake FetchContent:
- [GKlib](https://github.com/KarypisLab/GKlib) - Core library providing data structures and utilities to METIS
- [METIS](https://github.com/KarypisLab/METIS) - For graph partitioning
- [Scotch](https://gitlab.inria.fr/scotch/scotch) - For graph partitioning and sparse matrix ordering  
- [Robin Hood Hashing](https://github.com/martinus/robin-hood-hashing) - High-performance hash tables and containers

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

# <enable_multithreading>: 1 to enable, 0 to disable
./bin/dag-test <input.dot> <num_partitions> <enable_multithreading>

# Example
./bin/dag-test ../test/example.dot 4 1

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
./bin/dag-test random_dag.dot 4 1

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

## Output Format

The program outputs:

- Partition assignments for each vertex
- Edge cut value

## Performance

This tool was compared to the original [implementation](https://github.com/GT-TDAlab/dagP).
As input for both partitioning tools, a random DAG was generated using the [`rand-dag.cpp`](test/rand-dag.cpp) tool found in the [`test`](test) directory.
This tool generates a random DAG with a given node count and a ratio of edges to nodes (in %).
The generator creates deeper rather than wider DAGs, to mimic the structure of DAGs found in machine learning workflows. The testing script that runs the 2 partitioning tools with
the different inputs can be found in the [`test`](test) directory.

Both tools were compared based on the resulting edge cut, the imbalance ratio and the execution time.
The imbalance ratio is measured as the absolute difference of two percentages.
The percentage of the heaviest partition to the total weight minus the ideal percentage. For example, 
if we request 4 partitions from a graph with total weight 200 and the heaviest is partition has weight 
70, the imbalance is (70/200 = 35%) - (100/4 = 25%) = 10%.

For each input, each tool ran with all available parameter combinations implemented in this tool (clustering/bisection/refinement algorithm). The original implementation offers more algorithms that were not mentioned in the paper and were omitted from this tool.
The original implementation also offers additional parameterization options (such as node traversal strategies like random, DFS, or BFS) 
that were not implemented in this tool. For fair comparison, the original implementation was configured to match this tool's fixed 
approach where possible - for example, using topological ordering for node traversal. 

The results are presented per metric. For each metric and each input (DAG + partitions requested), the best performing parameter 
combination was considered for each tool.

The results can be found in the [`test`](test) directory, both comparative and absolute.

Finally, memory profiling was performed for this implementation, using the Valgrind Massif tool. Results show memory footprint across different graph sizes, types and partition counts. Similarly to the other metrics, the best performing algorithm combination was used for the graphs.

**Note on reliability:** The original tool sometimes failed to produce a solution and crashed with varying error messages. This appears in the heatmaps as red boxes indicating crashes across all parameter combinations for a specific DAG and partition count. Additional crashes occurred but aren't visible in the heatmap when at least one parameter combination succeeded. This implementation remained stable throughout all experiments with one exception. Due to heavy multithreading, valgrind rarely failed to run and crashed. If you encounter crashes with this tool, please open an issue.

**Note on performance:** This tool supports multi-threaded execution, unlike the original implementation. All benchmark results were obtained with this tool running multi-threaded (utilizing all available CPU threads) while the original tool ran single-threaded by necessity. This gives this implementation a significant advantage in execution time comparisons. The multi-threading behavior is user-controllable (see `test.cpp` and the `RecursiveBisectioner` constructor). For problems where this tool doesn't spawn additional threads (typically large graphs with few partitions), single-threaded performance can be observed in the time difference plots, showing this tool would perform worse in pure single-threaded comparisons. To compare single-threaded performance, modify [line 606](test/run-compare.py#L606) in [`run-compare.py`](test/run-compare.py) to use "0" instead of "1".

### Reproducing Results

To reproduce the performance comparison results:

1. Build and install this project (see Building section).
2. Build the original [dagP implementation](https://github.com/GT-TDAlab/dagP). Use this [config.py](https://gist.github.com/pxanthopoulos/b7891ce34dbefda2ad3499470e35b6fc) for building dagP, as well as the same external libraries that were built by this project.
3. Compile the [driver](https://gist.github.com/pxanthopoulos/da18d9609d12eaa7b4b9923c962892e8) for dagP and link with the dagP library.
4. Run the comparison script `python3 test/run-compare.py --baseline-executable /path/to/dagp`.

To reproduce the memory profiling results:

1. Build and install this project (see Building section).
2. Run the comparison script `python3 test/run-memprof.py`.

Both scripts support additional options - use `-h` or `--help` to see all available options.

*Note: The comparison results and analysis tools are available in the [`test`](test) directory.*

## Experimental Setup

All experiments and benchmarks were conducted on the following system:

- **CPU**: Intel Core i7-11370H @ 3.30GHz (4 cores, 8 threads, max 4.8GHz)
- **Memory**: 16GB RAM 
- **Operating System**: Ubuntu 22.04.5 LTS (Jammy Jellyfish)
- **Kernel**: Linux 6.8.0-65-generic
- **Architecture**: x86_64

The system includes support for AVX-512 instruction sets and Intel VT-x virtualization technology. All performance measurements were taken on this configuration.

## License

CC BY-NC 4.0

## References

Based on the [paper](https://epubs.siam.org/doi/abs/10.1137/18M1176865):
"Multilevel Algorithms for Acyclic Partitioning of Directed Acyclic Graphs" by Herrmann et al.
