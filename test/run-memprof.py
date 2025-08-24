# Run current implementation for memory profiling
try:
    with (
        open(args.current_trace_file, "r+b") as f,
        open(args.subprocess_log, "a") as log_f,
        open(massif_txt_path, "w+") as massif_f,
    ):
        subprocess.run(
            [
                "timeout",
                "--foreground",
                args.timeout,
                "valgrind",
                "--tool=massif",
                "--stacks=yes",
                f"--massif-out-file={massif_out_path}",
                test_executable_path,
                str(p),
                dag_dot_path,
                algorithm.split("-")[0],
                algorithm.split("-")[1],
                algorithm.split("-")[2],
                "1",
            ],
            stdout=subprocess.DEVNULL,
            stderr=log_f,
            check=True,
        )
        subprocess.run(
            [
                "ms_print",
                massif_out_path,
            ],
            stdout=massif_f,
            stderr=log_f,
            check=True,
        )

        footprint = 0.0

        massif_f.seek(0)
        lines = massif_f.read().strip().split("\n")
        empty_count = 0
        for i, line in enumerate(lines):
            if line.strip() == "":
                empty_count += 1
                if empty_count == 2:
                    unit = lines[i + 1].strip().split()[-1]
                    peak = float(lines[i + 2].split("^")[0].strip())
                    multipliers = {
                        "B": 1 / 1024,
                        "KB": 1,
                        "MB": 1024,
                        "GB": 1024 * 1024,
                    }
                    footprint = peak * multipliers[unit]
                    break
        f.seek(-1, 2)
        f.write(f",{footprint}\n".encode())
        f.flush()
        print(footprint)
except subprocess.CalledProcessError as e:
    with open(args.current_trace_file, "r+b") as f:
        f.seek(-1, 2)
        f.write(",\n".encode())
        f.flush()
    ret = e.returncode
    current_failures += 1
    error_msg = f"Error: current implementation failed for size {size}, ratio {ratio}, partitions {p} on iteration {run}"
    timeout_msg = "TIMED OUT" if ret == 124 else ""
    with open(args.subprocess_log, "a") as log_f:
        log_f.write(f"{error_msg}\n")
        if timeout_msg:
            log_f.write(f"{timeout_msg}\n")

# Run baseline implementation for memory profiling
try:
    with (
        open(args.baseline_trace_file, "r+b") as f,
        open(args.subprocess_log, "a") as log_f,
        open(massif_txt_path, "w+") as massif_f,
    ):
        subprocess.run(
            [
                "timeout",
                "--foreground",
                args.timeout,
                "valgrind",
                "--tool=massif",
                "--stacks=yes",
                f"--massif-out-file={massif_out_path}",
                args.baseline_executable,
                dag_dot_path,
                str(p),
                DAGP_CLUSTERING_MAP[algorithm.split("-")[0]],
                DAGP_BISECTION_MAP[algorithm.split("-")[1]],
                DAGP_REFINEMENT_MAP[algorithm.split("-")[2]],
            ],
            stdout=subprocess.DEVNULL,
            stderr=log_f,
            check=True,
        )
        subprocess.run(
            [
                "ms_print",
                massif_out_path,
            ],
            stdout=massif_f,
            stderr=log_f,
            check=True,
        )

        footprint = 0.0

        massif_f.seek(0)
        lines = massif_f.read().strip().split("\n")
        empty_count = 0
        for i, line in enumerate(lines):
            if line.strip() == "":
                empty_count += 1
                if empty_count == 2:
                    unit = lines[i + 1].strip().split()[-1]
                    peak = float(lines[i + 2].split("^")[0].strip())
                    multipliers = {
                        "B": 1 / 1024,
                        "KB": 1,
                        "MB": 1024,
                        "GB": 1024 * 1024,
                    }
                    footprint = peak * multipliers[unit]
                    break
        f.seek(-1, 2)
        f.write(f",{footprint}\n".encode())
        f.flush()
        print(footprint)
except subprocess.CalledProcessError as e:
    with open(args.baseline_trace_file, "r+b") as f:
        f.seek(-1, 2)
        f.write(",\n".encode())
        f.flush()
    ret = e.returncode
    current_failures += 1
    error_msg = f"Error: baseline implementation failed for size {size}, ratio {ratio}, partitions {p} on iteration {run}"
    timeout_msg = "TIMED OUT" if ret == 124 else ""
    with open(args.subprocess_log, "a") as log_f:
        log_f.write(f"{error_msg}\n")
        if timeout_msg:
            log_f.write(f"{timeout_msg}\n")
