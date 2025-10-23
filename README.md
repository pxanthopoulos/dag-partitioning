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

## Requirements

- A C++17 compatible compiler.
- CMake 3.28 or later.
- (Optional) Valgrind for memory profiling.

## Building the Project

This project uses CMake to automatically download and build its dependencies (GKlib, METIS, Scotch, and Robin Hood Hashing).

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

Based on the [paper](https://epubs.siam.org/doi/abs/10.1137/18M1176865):
"Multilevel Algorithms for Acyclic Partitioning of Directed Acyclic Graphs" by Herrmann et al.
