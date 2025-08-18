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

After building, you can use the `test` executable to partition DAGs:

```bash
# From the install directory
./bin/test input.dot num_partitions 1

# Example
./bin/test ../test/example.dot 4 1

# For detailed usage
./bin/test
```

### Generating Random DAGs

Use the `rand-dag` tool to generate test graphs:

```bash
# From the install directory
# Generate a random DAG with 100 nodes and 150% edge density
./bin/rand-dag 1000 150 0

# This creates a DOT format graph that can be used with the test executable
./bin/rand-dag 1000 150 0 random_dag.dot
./bin/test random_dag.dot 4 1

# For detailed usage
./bin/rand-dag
```

### Library Usage

For programmatic usage, check the test executable source in `[test/test.cpp](test/test.cpp)`.

## Input Format

The input DAG should be provided in the dot format as follows:

- First line: comment of the size N of the graph (number of vertices)
- Following N lines for vertices: *vertexname*[weight=*vertexweight*];
- Following lines for edges: *from*->*to*[weight=*edgeweight*];

Example:

```
// size=11
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
As input for the 2 partitioning tools, a random DAG was generated using the `rand-dag.cpp` tool found in the `test` directory.
This tool generates a random DAG with a given node count and a ratio of edges to nodes (in %).
The generator favours length over width, to mimic DAGs found in ML. The testing script that runs the 2 partitioning tools with
the different inputs can be found in the `test` directory.

The 2 tools were compared based on the resulting edge cut, the imbalance ratio and the execution time. 
The imbalance ratio is measured as the absolute difference of two percentages.
The percentage of the heaviest partition to the total weight minus the ideal percentage. For example, 
if we request 4 partitions from a graph with total weight 200 and the heaviest is partition has weight 
70, the imbalance is (70/200 = 35%) - (100/4 = 25%) = 10%.

For each input, each tool ran with all available parameter combinations (clustering/bisection/refinement algorithm). 
The original implementation offers additional parameterization options (such as node traversal strategies like random, DFS, or BFS) 
that were not implemented in this tool. For fair comparison, the original implementation was configured to match this tool's fixed 
approach where possible - for example, using topological ordering for node traversal. For each input, the best performing parameter 
combination was used for each tool.

The results can be found in the `test` directory.

Note on reliability: Sometimes, the original tool failed to produce a solution and crashed, with varying error messages. This can be 
seen in the heatmaps as white boxes (not including those at the bottom and left) indicating the original tool 
crashed across all parameter combinations. Additional crashes occurred but aren't apparent in the heatmap when at least one parameter 
combination produced a solution. This implementation remained stable throughout all experiments. If you
encounter any cases where this tool crashes, please feel free to open an issue.

## License

CC BY-NC 4.0

## References

Based on the [paper](https://epubs.siam.org/doi/abs/10.1137/18M1176865):
"Multilevel Algorithms for Acyclic Partitioning of Directed Acyclic Graphs" by Herrmann et al.
