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