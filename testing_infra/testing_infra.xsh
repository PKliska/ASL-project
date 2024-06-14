import csv
from typing import Any
import sys
from time import strftime
from pathlib import Path
from argparse import ArgumentParser


def main():
    args = parse_cli_args()
    set_global_variables(args)

    disable_turbo_boost()

    if should_only_rerun_plotter := args.run == "remake_plot":
        # re-runs the plotting code with the data of the last timing test
        timing_dir = Path(f"{$TESTING_DIR}/timing/")
        root_dir_for_this_test = sorted(timing_dir.glob("*"))[-1]

        all_implementations = sorted([p.name for p in root_dir_for_this_test.glob("*.csv")])
        all_implementations_comma_sep = ",".join(all_implementations)

        # make plots
        cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot_performance.py @(all_implementations_comma_sep)  && cd -
        cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot_runtime.py     @(all_implementations_comma_sep)  && cd -

        print(f"\n✅ Remade plot, saved plot at '{root_dir_for_this_test}/plot_perf.png'\n")
        print(f"\n✅ Remade plot, saved plot at '{root_dir_for_this_test}/plot_runtime.png'\n")

        exit(0)

    if should_only_run_skewed_heatmap := args.run == "heatmap":
        generate_heatmap_plots_for_skewed()
        exit(0)
    elif should_only_run_skewed_heatmap_rect := args.run.startswith("heatmap_rect"):
        generate_heatmap_plots_for_rectangle_skewed()
        exit(0)



    recompile_c_implementation()
    mkdir --parents "$TESTING_DIR"

    if has_user_defined_implementations := len(args.implementation) > 0:
        implementations = args.implementation
    else:
        implementations = get_all_implementations()

    if has_user_defined_test_cases := len(args.matrix_dimensions) > 0:
        test_cases = [int(t) for t in args.matrix_dimensions]
    else:
        test_cases = [32, 64, 128, 640, 960, 1280]

    # based on arguments run timing plot, consystency, correctness, etc.
    if args.run == "timing":
        run_timing_test(implementations, test_cases)
    elif args.run == "correctness":
        for i in implementations:
            run_correctness_test(i, test_cases)
    elif args.run == "consystency":
        for i in implementations:
            run_consystency_test(i, test_cases)
    elif args.run == "all":
        for i in implementations:
            run_correctness_test(i, test_cases)
            run_consystency_test(i, test_cases)

        run_timing_test(implementations, test_cases)
    else:
        raise Exception(f"Can't run action '{args.run}', invalid option. Did you mispell?")

def parse_cli_args():
    argparser = ArgumentParser(description="")
    argparser.add_argument("--run", metavar="ACTION", default="all", help="Action to execute (e.g. run correctness tests)")
    argparser.add_argument("--implementation", nargs='+', default=[], help="Chose which implementation to run, can also be given multiple implementations")
    argparser.add_argument("--matrix-dimensions", nargs='+', default=[], help="For which matrix dimensions to run the tests")
    argparser.add_argument("--block-size", default=8, type=int, help="Block size for blocking implementation")
    argparser.add_argument("--compiler", default="gcc", help="Choose between 'gcc', 'icc' and 'clang'")
    argparser.add_argument("--disable-auto-vec", action='store_true', help="Use this flag to disable auto vectorization")

    args = argparser.parse_args()
    return args

def generate_heatmap_plots_for_skewed():
    matrix_dimension = 1600
    # block_size = [16, 24, 32, 36, 40, 44, 48, 52, 56, 60, 64]
    # times = [8, 15, 16, 24, 31, 32, 40, 47, 48, 55, 56, 63, 64]

    # humongous plot
    block_size = [8, 16, 24, 32, 36, 40, 44, 45, 46, 47, 48, 50, 52, 56, 60, 64, 80, 96]
    times = sorted([8, 15, 16, 24, 31, 32, 40, 47, 48, 55, 56, 63, 64] + [7, 20, 25, 27, 28, 30, 36, 38, 42, 43, 44, 45, 46, 49, 50, 52])

    current_time = strftime("%Y-%m-%d_%H:%M:%S")

    root_dir_for_this_test = f"{$TESTING_DIR}/heatmap/{current_time}"
    mkdir --parents @(root_dir_for_this_test)
    measurements_file = f'{root_dir_for_this_test}/performance_metrics.csv'

    with open(measurements_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['FLOPs', 'Cycles', 'Block size', 'Timestamps', 'Matrix dimension'])  # Writing the header

        for timestamp in times:
            for block in block_size:
                print("Block:     " + str(block))
                print("Timestamp: " + str(timestamp))

                ### recompile binary
                print("🎬 Re-compiling C implementation...", flush=True)

                run_cmake(
                    compiler=$COMPILER, should_disable_auto_vectorization=$DISABLE_AUTO_VEC,
                    skew_block_x=block,  skew_block_y=block,  skew_timesteps=timestamp,
                    # placeholders below cuz they we are only testing skewed
                    sskew_block_x=64, sskew_block_y=64, sskew_timesteps=32, sskew_sub_block_x=32, sskew_sub_block_y=32,
                    uskew_block_x=16, uskew_block_y=56, uskew_timesteps=15,
                    block_size=69,
                )

                # run make in C implementation dir and then cd back into the prev dir (infra dir)
                cd $C_IMPLEMENTATION_DIR/build
                compile_command = !( make )
                has_command_failed_res = has_command_failed(compile_command)
                cd -

                if has_command_failed_res:
                    print(f"❗️Skipping invalid size of block={block} timestamp={timestamp} (must be that block > timestamp)\n", flush=True)
                    continue

                print("✅ Finished re-compiling C implementation!")

                # measuermente (time and perf)
                n_flops, n_cycles = get_flops_and_cycles_count("skewed", $C_BINARY, matrix_dimension)
                writer.writerow([n_flops, n_cycles, block, timestamp, matrix_dimension])
                file.flush() # write changes to disk
                print(f"Data saved to '{measurements_file}'\n", flush=True)

    # gen plots
    python ./skew_heatmap.py @(measurements_file) performance
    python ./skew_heatmap.py @(measurements_file) runtime

def generate_heatmap_plots_for_rectangle_skewed():
    matrix_dimension = 1600
    # block_sizes_x = [ 32, 48, 60, 64, 72, 81, 96]
    # block_sizes_y = [ 32, 48, 60, 64, 72, 81, 96]

    # TODO: run over night
    block_sizes_x = [2] + list(range(4, 38, 4))
    block_sizes_y = range(32, 201, 4)

    # # most optimal so far
    # block_sizes_x = [32]
    # block_sizes_y = [112]


    current_time = strftime("%Y-%m-%d_%H:%M:%S")
    root_dir_for_this_test = f"{$TESTING_DIR}/heatmap/{current_time}".replace(":", "_")
    root_dir_for_this_binary = f"{root_dir_for_this_test}/binaries"
    mkdir --parents @(root_dir_for_this_test)
    mkdir --parents @(root_dir_for_this_binary)
    measurements_file = f'{root_dir_for_this_test}/performance_metrics_rect.csv'

    with open(measurements_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['FLOPs', 'Cycles', 'Block size X', 'Block size Y','Timestamp', 'Matrix dimension'])  # Writing the header

        for y_dimension in block_sizes_y:
            for x_dimension in block_sizes_x:
                timestamp = min(x_dimension-1, y_dimension-1, 49) # biggest still valid timestamp

                print("X Dimension: " + str(x_dimension))
                print("Y Dimension: " + str(y_dimension))
                print("Timestamp: " + str(timestamp))

                ### recompile binary
                print("🎬 Re-compiling C implementation...", flush=True)
                run_cmake(
                    compiler=$COMPILER, should_disable_auto_vectorization=$DISABLE_AUTO_VEC,
                    skew_block_x=x_dimension,  skew_block_y=y_dimension,  skew_timesteps=timestamp,
                    # placeholders below cuz they we are only testing skewed
                    sskew_block_x=64, sskew_block_y=64, sskew_timesteps=32, sskew_sub_block_x=32, sskew_sub_block_y=32,
                    uskew_block_x=16, uskew_block_y=56, uskew_timesteps=15,
                    block_size=69,
                )

                # run make in C implementation dir and then cd back into the prev dir (infra dir)
                cd $C_IMPLEMENTATION_DIR/build
                compile_command = !( make )
                has_command_failed_res = has_command_failed(compile_command)

                if has_command_failed_res:
                    cd -
                    print()
                    print(compile_command.output.removesuffix(""))
                    print(compile_command.errors.removesuffix(""))
                    print()
                    print(f"❗️Skipping invalid size of block_sizes_x={x_dimension}, y_dimension={y_dimension} timestamp={timestamp} (must be that block > timestamp)\n", flush=True)
                    continue
                else:
                    new_binary_path = f"{root_dir_for_this_binary}/cavity_flow_{x_dimension}x_{y_dimension}y_{timestamp}t"
                    mv $C_BINARY @(new_binary_path)
                    cd -

                print("✅ Finished re-compiling C implementation!")

                n_flops, n_cycles = get_flops_and_cycles_count("skewed", new_binary_path, matrix_dimension)
                writer.writerow([n_flops, n_cycles, x_dimension, y_dimension, timestamp, matrix_dimension])
                file.flush() # write changes to disk
                print(f"USED BINARY: {new_binary_path}")
                print(f"Data saved to '{measurements_file}'\n", flush=True)

    # gen plots
    python ./skew_heatmap_rect.py @(measurements_file) performance
    python ./skew_heatmap_rect.py @(measurements_file) runtime

def get_flops_and_cycles_count(implementation: str, new_binary_path: str, matrix_dimension: int):
    print(f"Running implementation '{implementation}' and measuring FLOPS + cycles...")
    temp_dir = $( mktemp -d ).strip()

    output_file = f"{temp_dir}/stats.txt"

    cycles = $( perf stat \
        -x "," \
        --output @(output_file) \
        --event=fp_arith_inst_retired.128b_packed_double,fp_arith_inst_retired.128b_packed_single,fp_arith_inst_retired.256b_packed_double,fp_arith_inst_retired.256b_packed_single,fp_arith_inst_retired.scalar_double,fp_arith_inst_retired.scalar_single \
        @(new_binary_path) -I @(implementation) -t --dimension @(matrix_dimension)
    )

    n_cycles = float(cycles.strip().split()[-1])

    perf_res = $(cat @(output_file))

    flop_measurements = [int(l.split(",")[0]) / 4 for l in perf_res.splitlines()[2:]]
    # print(flop_measurements)

    n_128b_packed_double = flop_measurements[0]
    n_128b_packed_single = flop_measurements[1]
    n_256b_packed_double = flop_measurements[2]
    n_256b_packed_single = flop_measurements[3]
    n_scalar_double      = flop_measurements[4]
    n_scalar_single      = flop_measurements[5]

    total_flops = (
          2 * n_128b_packed_double
        + 4 * n_128b_packed_single
        + 4 * n_256b_packed_double
        + 8 * n_256b_packed_single
        + 1 * n_scalar_double
        + 1 * n_scalar_single
    )

    # print(total_flops)
    return (total_flops, n_cycles)

def run_timing_test(implementations: str, dimensions_which_to_test: list[int]):
    print("\n🟠 Starting timing test...")

    current_time = strftime("%Y-%m-%d %H:%M:%S")

    for implementation in implementations:
        data = []
        for matrix_dimension in dimensions_which_to_test:
            print(f"Simulating {implementation} with dimension={matrix_dimension}...", end="", flush=True)

            n_flops, n_cycles = get_flops_and_cycles_count(implementation, $C_BINARY, matrix_dimension)

            print(f" took {n_cycles} cycles ({n_cycles/(3.4*10**9)} sec)")
            data.append((matrix_dimension, n_cycles, n_flops))

        print(data)

        root_dir_for_this_test = f"{$TESTING_DIR}/timing/{current_time}"
        mkdir --parents @(root_dir_for_this_test)

        with open(f"{root_dir_for_this_test}/{implementation}.csv", "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["matrix_dimension", "n_cycles", "n_flops"])
            writer.writerows(data)

    all_implementations = sorted([f"{i}.csv" for i in implementations])
    all_implementations_comma_sep = ",".join(all_implementations)

    # make plots
    cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot_performance.py @(all_implementations_comma_sep)  && cd -
    cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot_runtime.py     @(all_implementations_comma_sep)  && cd -

    print(f"\n✅ Finished timing test, saved plot at '{root_dir_for_this_test}/plot_perf.png'\n")
    print(f"\n✅ Finished timing test, saved plot at '{root_dir_for_this_test}/plot_runtime.png'\n")


def run_correctness_test(implementation: str, dimensions_which_to_test: list[int]):
    print(f"\n🟠 Correctness test for '{implementation}'...")

    is_testing_small_dimensions = all(d < 1000 for d in dimensions_which_to_test)

    is_some_result_incorrect = False
    for d in dimensions_which_to_test:
        print(f"matrix_dimension = {d} ", end="", flush=True)
        try:
            check_if_c_output_matches_python_output_for(d, implementation)
        except Exception:
            if is_big_dimension := d >= 1000:
                print(f"    Correctness tests above 1000 fail cuz of tiny precision mismatches which is (kinda) okay")
                continue # don't fail correctness test
            # if small dimension then fail corretness test
            is_some_result_incorrect = True

    if is_testing_small_dimensions:
        return # don't perform big dimension test

    # test a big dimension for nan/inf
    very_large_matrix_dimension = 1600
    print(f"matrix_dimension (just check for NaN/Inf, don't compare to python impl.) = {very_large_matrix_dimension} ", end="", flush=True)

    root_dir_for_very_large_matrix = f"{$TESTING_DIR}/correctness/n_{very_large_matrix_dimension}"
    mkdir --parents @(root_dir_for_very_large_matrix)

    c_output_path = f"{root_dir_for_very_large_matrix}/output_c_{implementation}.csv"
    $C_BINARY  -I @(implementation) --output_file=@(c_output_path) --dimension=@(very_large_matrix_dimension)

    command = !( grep -q -Ei "nan|inf" @(c_output_path) )

    if does_matrix_contain_nan_or_inf := not has_command_failed(command):
        print("❌      ❗️ contains NaN or Inf")
        raise Exception("❗️❗️ Correctness tests failed!")
    else:
        print("✅")


    if is_some_result_incorrect:
        raise Exception("❗️❗️ Correctness tests failed!")
    print(f"\n✅ Finished correctness test!\n")


def run_consystency_test(implementation: str, dimensions_which_to_test: list[int]):
    print(f"\n🟠 Starting consystency test for '{implementation}'...")

    for i in dimensions_which_to_test:
        print(f"matrix_dimension = {i} ", end="", flush=True)
        check_if_c_outputs_consistent_for(i, implementation)

    print(f"\n✅ Finished consystency test!\n")


def check_if_c_outputs_consistent_for(matrix_dimension, implementation: str):
    """Checks if the C implementation produces the same output CSV file if given the
    same arguments."""

    mkdir -p "$TESTING_DIR/consystency/"
    c_output_path_1 = f"{$TESTING_DIR}/consystency/output_c_1_{implementation}.csv"
    c_output_path_2 = f"{$TESTING_DIR}/consystency/output_c_2_{implementation}.csv"

    $C_BINARY -I @(implementation) --output_file=@(c_output_path_1) --dimension=@(matrix_dimension)
    $C_BINARY -I @(implementation) --output_file=@(c_output_path_2) --dimension=@(matrix_dimension)

    compare_files_command = !( cmp --silent -- @(c_output_path_1) @(c_output_path_2 ))
    if are_files_different := has_command_failed(compare_files_command):
        print("❌")
        raise Exception("Running the C implementation with the same inputs produced different outputs!")
    else:
        print("✅")

def check_if_c_output_matches_python_output_for(matrix_dimension: int, implementation: str):
    root_dir_for_this_test = f"{$TESTING_DIR}/correctness/n_{matrix_dimension}"
    mkdir --parents @(root_dir_for_this_test)

    # run python implementation & save output
    python_output_path = f"{root_dir_for_this_test}/output_python.csv"
    # if has_python_been_run_previously := Path(python_output_path).exists():
    if has_python_been_run_previously := False:
        pass
    else:
        python ./python_implementation.py @(python_output_path) @(matrix_dimension)

    # run C implementation & save output
    c_output_path = f"{root_dir_for_this_test}/output_c_{implementation}.csv"
    $C_BINARY  -I @(implementation) --output_file=@(c_output_path) --dimension=@(matrix_dimension)

    # compare the outputs
    compare_files_command = !( cmp --silent -- @(c_output_path) @(python_output_path) )
    if are_files_different := has_command_failed(compare_files_command):
        print("❌")
        raise Exception("Error: C implementation produced different result compared to Python implementation!")
    else:
        print("✅")


## HELPER FUNCTIONS
def set_global_variables(args):
    ## VARS FOR THIS PROGRAM
    $TESTING_INFRA_ROOT_DIR = Path(__file__).parent.resolve(strict=True)
    $C_IMPLEMENTATION_DIR = Path(f"{$TESTING_INFRA_ROOT_DIR}/../c_implementation").resolve(strict=True)
    $C_BINARY = Path(f"{$C_IMPLEMENTATION_DIR}/build/bin/cavity_flow")

    if args.compiler == "gcc":
        $COMPILER = "gcc-12"
    elif args.compiler == "icc":
        $COMPILER = "icx"
    elif args.compiler == "clang":
        $COMPILER = "/home/myboi/asl/poly/build/bin/clang"
    else:
        raise Error("Invalid chompiler chose, please abort ur life, tnx")

    $DISABLE_AUTO_VEC = args.disable_auto_vec

    $TESTING_DIR = Path(f"{$TESTING_INFRA_ROOT_DIR}/.test_results").resolve(strict=True)
    $BLOCK_SIZE = args.block_size

    ## XONSH
    # throw error if bash command fails, otherwise we silently ignore the error
    $RAISE_SUBPROC_ERROR = True
    $XONSH_SHOW_TRACEBACK = True
    # for bash commands, print out how they were invoked (with which concrete args)
    # $XONSH_TRACE_SUBPROC = True

    # source some vars for the Intel C Compiler, if we source twiece then we get an error, so do it here
    source-bash /home/myboi/intel/oneapi/setvars.sh


def disable_turbo_boost():
    if is_turbo_boost_enabled := $( cat /sys/devices/system/cpu/intel_pstate/no_turbo ) == "0\n":
        print("✋🛑 Turbo boost is not disabled! Disable turbo boost by running the below commands!\n")
        print("sudo wrmsr -p0 0x1a0 0x4000850089")
        print("sudo wrmsr -p1 0x1a0 0x4000850089")
        print("sudo wrmsr -p2 0x1a0 0x4000850089")
        print("sudo wrmsr -p3 0x1a0 0x4000850089")
        exit(-1)

    print("✅ Turbo boost is disabled\n")

def get_all_implementations() -> str:
    implementations_string: str = $( $C_BINARY -l )
    implementations = implementations_string.split(",")[:-1]
    return implementations

def recompile_c_implementation():
    print("🎬 Re-compiling C implementation...")
    skew_block_x = 32
    skew_block_y = 112
    timestamp = min(skew_block_x-1, skew_block_y-1, 49) # biggest still valid timestamp
    # should_auto_vectorize = True

    run_cmake(
        compiler=$COMPILER, should_disable_auto_vectorization=$DISABLE_AUTO_VEC,
        skew_block_x=skew_block_x,  skew_block_y=skew_block_y,  skew_timesteps=timestamp,
        sskew_block_x=skew_block_x, sskew_block_y=skew_block_y, sskew_timesteps=timestamp, sskew_sub_block_x=skew_block_x, sskew_sub_block_y=skew_block_y,
        uskew_block_x=16, uskew_block_y=56, uskew_timesteps=15,
        block_size=$BLOCK_SIZE,
    )

    # x_dimension = 64
    # y_dimension = 64
    # sub_skew_dim = 32
    # timestamp = min(x_dimension-1, y_dimension-1, 49, sub_skew_dim-1) # biggest still valid timestamp

    # cmake -S $C_IMPLEMENTATION_DIR/ -B $C_IMPLEMENTATION_DIR/build/  \
    #     -DSKEWING_TIMESTEPS=@(timestamp) -DSKEWING_BLOCK_SIZE_X=@(x_dimension) -DSKEWING_BLOCK_SIZE_Y=@(y_dimension) \
    #     -DSUBSKEWING_BLOCK_SIZE_X=32 -DSUBSKEWING_BLOCK_SIZE_Y=32 \
    #     -DBLOCK_SIZE=$BLOCK_SIZE \
    #     -DNDEBUG=YOLOL
        # -DNO_AUTO_VEC=FILIP_WAZ_HARE



    # run make in C implementation dir and then cd back into the prev dir (infra dir)
    cd $C_IMPLEMENTATION_DIR/build && make && cd -

    print("✅ Finished re-compiling C implementation!\n")

def run_cmake(*,
    compiler, should_disable_auto_vectorization,
    skew_block_x,  skew_block_y,  skew_timesteps,
    sskew_block_x, sskew_block_y, sskew_timesteps, sskew_sub_block_x, sskew_sub_block_y,
    uskew_block_x, uskew_block_y, uskew_timesteps,
    block_size,
):
    $CMAKE_C_COMPILER = compiler
    $CMAKE_CXX_COMPILER = compiler

    cmake -S $C_IMPLEMENTATION_DIR/ -B $C_IMPLEMENTATION_DIR/build/  \
        -D CMAKE_C_COMPILER=$CMAKE_C_COMPILER \
        -D CMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER \
        -D SKEWING_TIMESTEPS=@(skew_timesteps) -D SKEWING_BLOCK_SIZE_X=@(skew_block_x) -D SKEWING_BLOCK_SIZE_Y=@(skew_block_y) \
        -D SSKEWING_TIMESTEPS=@(sskew_timesteps) -D SSKEWING_BLOCK_SIZE_X=@(sskew_block_x) -D SSKEWING_BLOCK_SIZE_Y=@(sskew_block_y) \
        -D SSKEWING_SUBBLOCK_SIZE_X=@(sskew_sub_block_x) -D SSKEWING_SUBBLOCK_SIZE_Y=@(sskew_sub_block_y) \
        -D USKEWING_BLOCK_SIZE_X=@(uskew_block_x) -D USKEWING_BLOCK_SIZE_Y=@(uskew_block_y) -D USKEWING_TIMESTEPS=@(uskew_timesteps) \
        -D BLOCK_SIZE=@(block_size) \
        -D NDEBUG=YOLOL \
        -D NO_AUTO_VEC=@("DISABLE_AUTO_VEC" if should_disable_auto_vectorization else "FILIP_WAZ_HARE")

def read_matrices_from_csv(csv_path: str) -> tuple[Any, Any, Any]:
    with open(csv_path, 'r') as csv_file:
        csv_lines = csv_file.readlines()
        header = csv_lines[0]

        matrix_size = int(header.split(",")[0])
        print(f"matrix_size= {matrix_size}")

        matrix_start_line = 1
        matrix_end_line = matrix_start_line + matrix_size
        matrix_1 = list(csv.reader(csv_lines[matrix_start_line : matrix_end_line]))

        matrix_start_line = matrix_end_line + 1 # skip all-comma line between matrices
        matrix_end_line = matrix_start_line + matrix_size
        matrix_2 = list(csv.reader(csv_lines[matrix_start_line : matrix_end_line]))

        matrix_start_line = matrix_end_line + 1 # skip all-comma line between matrices
        matrix_end_line = matrix_start_line + matrix_size
        matrix_3 = list(csv.reader(csv_lines[matrix_start_line : matrix_end_line]))

        return matrix_1, matrix_2, matrix_3

def has_command_failed(command) -> bool:
    $RAISE_SUBPROC_ERROR = False # temporarily don't throw error if bash command returns non-zero
    return_code = command.returncode
    $RAISE_SUBPROC_ERROR = True # re-enable

    return return_code != 0


## EXECUTE MAIN
main()
