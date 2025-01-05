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
- CMake 3.29 or later

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

```cpp
// Example usage code will be added here
```

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
- Partition sizes

## License

CC BY-NC 4.0

## References

Based on the [paper](https://epubs.siam.org/doi/abs/10.1137/18M1176865):
"Multilevel Algorithms for Acyclic Partitioning of Directed Acyclic Graphs" by Herrmann et al.

## Todo

- [x] Implement coarsening/clustering phase
- [x] Implement initial partitioning phase
- [ ] Implement refinement phase
- [ ] Add comprehensive tests
- [ ] Add benchmarking suite
- [ ] Complete documentation
- [ ] Add examples
