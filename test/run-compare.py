import subprocess
import argparse
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math

ALGORITHMS = [
    "FORB-UNDIRMETIS-BOUNDARYFM",
    "FORB-UNDIRMETIS-BOUNDARYKL",
    "FORB-UNDIRMETIS-MIXED",
    "FORB-GGG-BOUNDARYFM",
    "FORB-GGG-BOUNDARYKL",
    "FORB-GGG-MIXED",
    "CYC-UNDIRMETIS-BOUNDARYFM",
    "CYC-UNDIRMETIS-BOUNDARYKL",
    "CYC-UNDIRMETIS-MIXED",
    "CYC-GGG-BOUNDARYFM",
    "CYC-GGG-BOUNDARYKL",
    "CYC-GGG-MIXED",
    "HYB-UNDIRMETIS-BOUNDARYFM",
    "HYB-UNDIRMETIS-BOUNDARYKL",
    "HYB-UNDIRMETIS-MIXED",
    "HYB-GGG-BOUNDARYFM",
    "HYB-GGG-BOUNDARYKL",
    "HYB-GGG-MIXED",
]

DAGP_CLUSTERING_MAP = {"FORB": "0", "CYC": "1", "HYB": "2"}
DAGP_BISECTION_MAP = {"UNDIRMETIS": "0", "GGG": "1"}
DAGP_REFINEMENT_MAP = {"BOUNDARYFM": "0", "BOUNDARYKL": "1", "MIXED": "2"}

metrics = ["edge_cut", "imbalance", "time"]
metric_dict = {
    "edge_cut": "Edge Cut Difference %",
    "imbalance": "Imbalance Ratio Difference",
    "time": "Execution Time Difference %",
}
titles_diff = {
    "edge_cut": "Edge Cut Difference Percentage (current Implementation - Original)",
    "imbalance": "Imbalance Ratio Difference (current Implementation - Original)",
    "time": "Execution Time Difference Percentage (current Implementation - Original)",
}
titles_absolute = {
    "edge_cut": "Edge Cut",
    "imbalance": "Imbalance Ratio",
    "time": "Execution Time",
}


def print_summary(
    iterations_run, total_planned, current_failures, baseline_failures, interrupted
):
    """Print test execution summary"""
    status = "INTERRUPTED" if interrupted else "COMPLETED"
    print(f"\n=== Test Summary ({status}) ===")
    print(f"Iterations run: {iterations_run} / {total_planned}")
    if interrupted:
        print(
            f"Progress: {(iterations_run / total_planned * 100):.1f}% complete when interrupted"
        )
    print(f"Current implementation failures: {current_failures}")
    print(f"Baseline implementation failures: {baseline_failures}")

    if iterations_run > 0:
        print(
            f"Success rate - Current: {(((iterations_run * 18) - current_failures) / (iterations_run * 18) * 100):.1f}%"
        )
        print(
            f"Success rate - Baseline: {(((iterations_run * 18) - baseline_failures) / (iterations_run * 18) * 100):.1f}%"
        )
    else:
        print("No iterations completed.")


def cleanup_files(args, dag_dot_path, node_mappings_current, node_mappings_dagp):
    """Clean up generated files"""
    cleanup_files_list = [dag_dot_path, node_mappings_current, node_mappings_dagp]
    if not args.keep_subprocess_log:
        cleanup_files_list.append(args.subprocess_log)

    for file_path in cleanup_files_list:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Cleaned up: {file_path}")
        except OSError as e:
            print(f"Warning: Could not remove {file_path}: {e}")


def validate_ratios(ratios_str):
    """Parse and validate ratios from comma-separated string"""
    try:
        ratios = [int(x.strip()) for x in ratios_str.split(",")]
        for ratio in ratios:
            if ratio < 100:
                raise argparse.ArgumentTypeError(f"Ratio {ratio} must be >= 100")
        return ratios
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Ratios must be comma-separated integers: {ratios_str}"
        )


def validate_sizes(sizes_str):
    """Parse and validate sizes from comma-separated string"""
    try:
        sizes = [int(x.strip()) for x in sizes_str.split(",")]
        for size in sizes:
            if size <= 0:
                raise argparse.ArgumentTypeError(
                    f"Size {size} must be a positive integer"
                )
        return sizes
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Sizes must be comma-separated positive integers: {sizes_str}"
        )


def validate_timeout(timeout_str):
    """Parse and validate timeout duration string"""
    import re

    # Match floating point number with optional suffix
    pattern = r"^(\d*\.?\d+)([smhd]?)$"
    match = re.match(pattern, timeout_str.strip())

    if not match:
        raise argparse.ArgumentTypeError(
            f"Invalid timeout format: {timeout_str}. Expected format: number[smhd]"
        )

    value, suffix = match.groups()
    try:
        value = float(value)
        if value < 0:
            raise argparse.ArgumentTypeError(
                f"Timeout value must be non-negative: {value}"
            )
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid timeout value: {value}")

    # Validate suffix (s=seconds, m=minutes, h=hours, d=days)
    if suffix and suffix not in "smhd":
        raise argparse.ArgumentTypeError(
            f"Invalid timeout suffix: {suffix}. Valid suffixes: s, m, h, d"
        )

    return timeout_str.strip()  # Return original format for timeout command


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare dag-partitioning implementations"
    )
    parser.add_argument(
        "--current-trace-file",
        default="./trace-current.txt",
        help="Output trace file for current implementation (default: ./trace-current.txt)",
    )
    parser.add_argument(
        "--baseline-trace-file",
        default="./trace-dagp.txt",
        help="Output trace file for baseline (dagp) implementation (default: ./trace-dagp.txt)",
    )
    parser.add_argument(
        "--dag-output-file",
        default="./dag.dot",
        help="Temporary DAG file generated by rand-dag (default: ./dag.dot)",
    )
    parser.add_argument(
        "--subprocess-log",
        default="./subprocesses.log",
        help="Log file for subprocess errors (default: ./subprocesses.log)",
    )
    parser.add_argument(
        "--keep-subprocess-log",
        action="store_true",
        help="Keep subprocess log file after execution (default: remove it)",
    )
    parser.add_argument(
        "--ratios",
        type=validate_ratios,
        default="100,150,200,250,300,500",
        help="Comma-separated ratios (integers >= 100) (default: 100,150,200,250,300,500)",
    )
    parser.add_argument(
        "--sizes",
        type=validate_sizes,
        default="20,50,100,150,200,500,1000,2000,5000,10000,50000",
        help="Comma-separated sizes (positive integers) (default: 20,50,100,150,200,500,1000,2000,5000,10000,50000)",
    )
    parser.add_argument(
        "--timeout",
        type=validate_timeout,
        default="5m",
        help="Timeout duration for each test (format: number[smhd], e.g., 5m, 30s, 1h) (default: 5m)",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=5,
        help="Number of runs per partition configuration (default: 5)",
    )
    parser.add_argument(
        "--max-base-partitions",
        type=int,
        default=5,
        help="Maximum value for base partition range (2 to this value inclusive). Power-of-2 sequence starts from next power of 2 after this limit (default: 5)",
    )
    parser.add_argument(
        "--baseline-executable",
        required=True,
        help="Path to baseline (dagp) executable for comparison",
    )
    parser.add_argument(
        "--plot-path", default="./", help="Output path for heatmaps (default: ./)"
    )
    return parser.parse_args()


def calculate_partitions(size, max_base):
    partitions = []

    # Add base partitions from 2 to max_base (inclusive), but don't exceed size
    for i in range(2, min(max_base + 1, size + 1)):
        partitions.append(i)

    # Find next power of 2 after max_base
    power = 1
    while power <= max_base:
        power *= 2

    # Add power-of-2 partitions starting from the next power of 2 after max_base
    sizediv2 = size // 2
    while power <= sizediv2:
        partitions.append(power)
        power *= 2

    return partitions


def calculate_total_iterations(ratios, sizes, runs, max_base):
    total = 0
    for ratio in ratios:
        for size in sizes:
            partitions = calculate_partitions(size, max_base)
            total += len(partitions) * runs
    return total


def generate_iteration_space(ratios, sizes, max_base):
    ratio_values = []
    size_values = []
    partitions_values = []
    algo_values = []

    for ratio in ratios:
        for size in sizes:
            partitions = calculate_partitions(size, max_base)
            for partition in partitions:
                for algo in ALGORITHMS:
                    ratio_values.append(ratio)
                    size_values.append(size)
                    partitions_values.append(partition)
                    algo_values.append(algo)

    return ratio_values, size_values, partitions_values, algo_values


def parse_trace_file(filepath, ratios, sizes, max_base, runs):
    with open(filepath, "r") as file:
        trace_data = file.read()

    data = {
        "run": [],
        "ratio": [],
        "size": [],
        "partitions": [],
        "algorithm": [],
        "edge_cut": [],
        "imbalance": [],
        "time": [],
    }

    current_run = None
    current_size = None
    current_ratio = None
    current_partitions = None

    for line in trace_data.strip().split("\n"):
        if line.startswith("Error"):
            continue
        elif line.startswith("Counter"):
            parts = line.split()
            current_run = int(parts[3].rstrip(","))
            current_size = int(parts[5].rstrip(","))
            current_ratio = int(parts[7].rstrip(","))
            current_partitions = int(parts[9].rstrip(","))
        else:
            clustering, bisection, refinement, metric1, metric2, metric3 = line.split(
                ","
            )
            data["run"].append(current_run)
            data["ratio"].append(current_ratio)
            data["size"].append(current_size)
            data["partitions"].append(current_partitions)
            data["algorithm"].append(f"{clustering}-{bisection}-{refinement}")
            data["edge_cut"].append(int(metric1) if metric1.strip() else np.nan)
            data["imbalance"].append(float(metric2) if metric2.strip() else np.nan)
            data["time"].append(int(metric3) if metric3.strip() else np.nan)

    df = pd.DataFrame(data)
    grouped = df.groupby(
        ["ratio", "size", "partitions", "algorithm"], as_index=False
    ).agg(
        {
            "edge_cut": "mean",
            "imbalance": "mean",
            "time": "mean",
        }
    )

    csv_filepath = os.path.splitext(filepath)[0] + ".csv"
    grouped.to_csv(csv_filepath, index=False)

    print(f"Processed {filepath}")
    print(
        f"  Total expected rows: {int(calculate_total_iterations(ratios, sizes, runs, max_base) * len(ALGORITHMS) / runs)}"
    )
    print(f"  Rows with data: {len(grouped.dropna())}")
    print(f"  Rows with NaN: {len(grouped) - len(grouped.dropna())}")
    print(f"  Output saved to: {csv_filepath}")
    print()


def generate_diff_df(df_current, df_dagp, metric):
    other_metrics = [m for m in metrics if m != metric]
    df_current_metric = df_current.drop(columns=other_metrics)
    df_dagp_metric = df_dagp.drop(columns=other_metrics)
    dfs = {"current": df_current_metric, "dagp": df_dagp_metric}
    for name, df in dfs.items():
        df["group"] = df.index // 18

        best_indices = []
        for _, group_df in df.groupby("group"):
            if group_df[metric].isna().all():
                best_indices.append(group_df.index[0])
            else:
                best_indices.append(group_df[metric].idxmin())

        dfs[name] = df.loc[best_indices]

    df_current_metric = dfs["current"].drop(columns=["algorithm", "group"])
    df_dagp_metric = dfs["dagp"].drop(columns=["algorithm", "group"])
    df_current_metric = df_current_metric.reset_index(drop=True)
    df_dagp_metric = df_dagp_metric.reset_index(drop=True)

    diff_df = df_current_metric.copy()
    if metric == "imbalance":
        diff_df[metric] = diff_df[metric] - df_dagp_metric[metric]
    else:
        diff_df[metric] = (
            (diff_df[metric] - df_dagp_metric[metric]) / df_dagp_metric[metric] * 100
        )

    return diff_df


def generate_absolute_df(df, metric):
    other_metrics = [m for m in metrics if m != metric]
    df_metric = df.drop(columns=other_metrics)
    df_metric["group"] = df_metric.index // 18

    best_indices = []
    for _, group_df in df_metric.groupby("group"):
        if group_df[metric].isna().all():
            best_indices.append(group_df.index[0])
        else:
            best_indices.append(group_df[metric].idxmin())

    df_metric = df_metric.loc[best_indices]

    df_metric = df_metric.drop(columns=["algorithm", "group"])
    df_metric = df_metric.reset_index(drop=True)

    return df_metric


def plot_df(diff_df, metric, plot_title, file_title, base_path, width, height):
    ratios = sorted(diff_df["ratio"].unique())

    fig, axes = plt.subplots(height, width, figsize=(20, 10))
    axes = axes.flatten()
    for idx, ratio in enumerate(ratios):
        data = diff_df[diff_df["ratio"] == ratio]
        pivot = data.pivot(index="size", columns="partitions", values=metric)
        pivot = pivot.iloc[::-1]

        sns.heatmap(
            pivot,
            ax=axes[idx],
            cmap="viridis",
            annot=False,
            cbar_kws={"label": metric_dict[metric]},
            square=True,
            linewidths=0.5,
            linecolor="black",
        )

        for i, size_val in enumerate(pivot.index):
            for j, part_val in enumerate(pivot.columns):
                combo_exists = (
                    len(
                        data[
                            (data["size"] == size_val)
                            & (data["partitions"] == part_val)
                        ]
                    )
                    > 0
                )
                has_nan_value = pd.isna(pivot.iloc[i, j])

                if combo_exists and has_nan_value:
                    axes[idx].add_patch(
                        plt.Rectangle(
                            (j, i),
                            1,
                            1,
                            fill=True,
                            facecolor="red",
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=10,
                        )
                    )

        axes[idx].set_xticklabels(axes[idx].get_xticklabels(), rotation=45, ha="right")
        axes[idx].set_title(f"Ratio: {ratio}")
        axes[idx].set_xlabel("Partitions")
        axes[idx].set_ylabel("Size")
    fig.suptitle(plot_title, fontsize=22)
    plt.tight_layout()
    plt.savefig(f"{base_path}/{file_title}", bbox_inches="tight", dpi=500)


def plot_all_diffs(current_trace_file, baseline_trace_file, plot_path, width, height):
    current_trace_csv = os.path.splitext(current_trace_file)[0] + ".csv"
    baseline_trace_csv = os.path.splitext(baseline_trace_file)[0] + ".csv"
    df_current = pd.read_csv(current_trace_csv)
    df_dagp = pd.read_csv(baseline_trace_csv)

    for metric in metrics:
        diff_df = generate_diff_df(df_current, df_dagp, metric)
        plot_df(
            diff_df,
            metric,
            titles_diff[metric],
            f"heatmap-diff-{metric}.png",
            plot_path,
            width,
            height,
        )
        df_current_metric = generate_absolute_df(df_current, metric)
        plot_df(
            df_current_metric,
            metric,
            titles_absolute[metric] + " (Current Implementation)",
            f"heatmap-current-{metric}.png",
            plot_path,
            width,
            height,
        )
        df_dagp_metric = generate_absolute_df(df_dagp, metric)
        plot_df(
            df_dagp_metric,
            metric,
            titles_absolute[metric] + " (Original Implementation)",
            f"heatmap-dagp-{metric}.png",
            plot_path,
            width,
            height,
        )


def main():
    args = parse_args()

    # # Get script directory and project root
    # script_dir = os.path.dirname(os.path.abspath(__file__))
    # # Go up one level from test/ to project root
    # project_root = os.path.dirname(script_dir)

    # # Paths relative to project root
    # rand_dag_path = os.path.join(project_root, "install", "bin", "rand-dag")
    # test_executable_path = os.path.join(project_root, "install", "bin", "test")
    # dag_dot_path = os.path.abspath(args.dag_output_file)

    # # Files created by executables that need cleanup
    # node_mappings_current = dag_dot_path + ".node-mappings.txt"
    # node_mappings_dagp = dag_dot_path + ".nodemappings"

    # # Clear trace files
    # open(args.current_trace_file, "w").close()
    # open(args.baseline_trace_file, "w").close()
    # open(args.subprocess_log, "w").close()

    # counter = 1
    # total_iterations = calculate_total_iterations(
    #     args.ratios, args.sizes, args.runs, args.max_base_partitions
    # )
    # current_failures = 0
    # baseline_failures = 0
    # iterations_run = 0
    # interrupted = False

    # # Create progress bar with colors
    # pbar = tqdm(
    #     total=total_iterations,
    #     desc="Tests",
    #     bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    #     colour="green",
    # )

    # try:
    #     for ratio in args.ratios:
    #         for size in args.sizes:
    #             partitions = calculate_partitions(size, args.max_base_partitions)
    #             for p in partitions:
    #                 # Run rand-dag
    #                 try:
    #                     with open(args.subprocess_log, "a") as log_f:
    #                         subprocess.run(
    #                             [
    #                                 rand_dag_path,
    #                                 str(size),
    #                                 str(ratio),
    #                                 "0",
    #                                 dag_dot_path,
    #                             ],
    #                             check=True,
    #                             stderr=log_f,
    #                         )
    #                 except subprocess.CalledProcessError:
    #                     print(
    #                         f"Error: rand-dag failed for size {size}, ratio {ratio} (check {args.subprocess_log})"
    #                     )
    #                     raise RuntimeError(
    #                         f"rand-dag failed for size {size}, ratio {ratio}"
    #                     )

    #                 for i in range(1, args.runs + 1):
    #                     pbar.set_description(f"C{counter} R{i} S{size} Ra{ratio} P{p}")
    #                     pbar.update(1)
    #                     iterations_run += 1
    #                     with open(args.current_trace_file, "a") as f:
    #                         f.write(
    #                             f"Counter {counter}, Run {i}, Size {size}, Ratio {ratio}, Partitions {p}\n"
    #                         )
    #                     with open(args.baseline_trace_file, "a") as f:
    #                         f.write(
    #                             f"Counter {counter}, Run {i}, Size {size}, Ratio {ratio}, Partitions {p}\n"
    #                         )
    #                     for algorithm in ALGORITHMS:
    #                         # Run current implementation with timeout
    #                         try:
    #                             with (
    #                                 open(args.current_trace_file, "a") as f,
    #                                 open(args.subprocess_log, "a") as log_f,
    #                             ):
    #                                 subprocess.run(
    #                                     [
    #                                         "timeout",
    #                                         "--foreground",
    #                                         args.timeout,
    #                                         test_executable_path,
    #                                         str(p),
    #                                         dag_dot_path,
    #                                         algorithm.split("-")[0],
    #                                         algorithm.split("-")[1],
    #                                         algorithm.split("-")[2],
    #                                         "1",
    #                                     ],
    #                                     stdout=f,
    #                                     stderr=log_f,
    #                                     check=True,
    #                                 )
    #                         except subprocess.CalledProcessError as e:
    #                             with open(args.current_trace_file, "a") as f:
    #                                 f.write(
    #                                     f"{algorithm.split('-')[0]},{algorithm.split('-')[1]},{algorithm.split('-')[2]},,,\n"
    #                                 )
    #                             ret = e.returncode
    #                             current_failures += 1
    #                             error_msg = f"Error: current implementation failed for size {size}, ratio {ratio}, partitions {p} on iteration {i}"
    #                             timeout_msg = "TIMED OUT" if ret == 124 else ""
    #                             with open(args.subprocess_log, "a") as log_f:
    #                                 log_f.write(f"{error_msg}\n")
    #                                 if timeout_msg:
    #                                     log_f.write(f"{timeout_msg}\n")

    #                         # Run baseline implementation with timeout
    #                         try:
    #                             with (
    #                                 open(args.baseline_trace_file, "a") as f,
    #                                 open(args.subprocess_log, "a") as log_f,
    #                             ):
    #                                 subprocess.run(
    #                                     [
    #                                         "timeout",
    #                                         "--foreground",
    #                                         args.timeout,
    #                                         args.baseline_executable,
    #                                         dag_dot_path,
    #                                         str(p),
    #                                         DAGP_CLUSTERING_MAP[
    #                                             algorithm.split("-")[0]
    #                                         ],
    #                                         DAGP_BISECTION_MAP[algorithm.split("-")[1]],
    #                                         DAGP_REFINEMENT_MAP[
    #                                             algorithm.split("-")[2]
    #                                         ],
    #                                     ],
    #                                     stdout=f,
    #                                     stderr=log_f,
    #                                     check=True,
    #                                 )
    #                         except subprocess.CalledProcessError as e:
    #                             with open(args.baseline_trace_file, "a") as f:
    #                                 f.write(
    #                                     f"{algorithm.split('-')[0]},{algorithm.split('-')[1]},{algorithm.split('-')[2]},,,\n"
    #                                 )
    #                             ret = e.returncode
    #                             baseline_failures += 1
    #                             error_msg = f"Error: baseline implementation failed for size {size}, ratio {ratio}, partitions {p} on iteration {i}"
    #                             timeout_msg = "TIMED OUT" if ret == 124 else ""
    #                             with open(args.subprocess_log, "a") as log_f:
    #                                 log_f.write(f"{error_msg}\n")
    #                                 if timeout_msg:
    #                                     log_f.write(f"{timeout_msg}\n")

    #                 counter += 1

    # except KeyboardInterrupt:
    #     interrupted = True
    #     print("\n\nInterrupted by user (Ctrl+C)")

    # finally:
    #     pbar.close()

    #     # Print summary and cleanup
    #     print_summary(
    #         iterations_run,
    #         total_iterations,
    #         current_failures,
    #         baseline_failures,
    #         interrupted,
    #     )
    #     cleanup_files(args, dag_dot_path, node_mappings_current, node_mappings_dagp)

    # Parse trace files to CSV
    print("Parsing trace files to CSV...")
    parse_trace_file(
        args.current_trace_file,
        args.ratios,
        args.sizes,
        args.max_base_partitions,
        args.runs,
    )
    parse_trace_file(
        args.baseline_trace_file,
        args.ratios,
        args.sizes,
        args.max_base_partitions,
        args.runs,
    )

    print("Generating difference plots...")
    plot_all_diffs(
        args.current_trace_file,
        args.baseline_trace_file,
        args.plot_path,
        math.ceil(len(args.ratios) / 2),
        2,
    )


if __name__ == "__main__":
    main()
