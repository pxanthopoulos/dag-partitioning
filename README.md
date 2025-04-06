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

## Requirements

- C++17 or later
- CMake 3.28 or later

## Building

This project requires the following dependencies which must be built separately before proceeding:
- [METIS](https://github.com/KarypisLab/METIS) - For graph partitioning
- [Scotch](https://gitlab.inria.fr/scotch/scotch) - For graph partitioning and sparse matrix ordering

### General Build Instructions

```bash
mkdir build
cd build
cmake -DMETIS_INCLUDE_DIR=/path/to/metis/include \
      -DSCOTCH_INCLUDE_DIR=/path/to/scotch/include \
      -DMETIS_LIBRARY=/path/to/metis/lib/libmetis.a \
      -DSCOTCH_LIBRARY=/path/to/scotch/lib/libscotch.a \
      -DSCOTCHERR_LIBRARY=/path/to/scotch/lib/libscotcherr.a \
      ..
make
```

### Example with Custom Paths

```bash
mkdir build
cd build
cmake -DMETIS_INCLUDE_DIR=/home/user/METIS/build/install/include \
      -DSCOTCH_INCLUDE_DIR=/home/user/scotch/build/src/include \
      -DMETIS_LIBRARY=/home/user/METIS/build/install/lib/libmetis.a \
      -DSCOTCH_LIBRARY=/home/user/scotch/build/lib/libscotch.a \
      -DSCOTCHERR_LIBRARY=/home/user/scotch/build/lib/libscotcherr.a \
      ..
make
```

## Usage

For example usage, see `main.cpp` at the `aio-executable` branch.

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
