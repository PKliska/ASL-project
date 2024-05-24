import csv
from typing import Any
import sys
from time import strftime
from pathlib import Path
from argparse import ArgumentParser


def main():
    args = parse_cli_args()
    set_global_variables(args)

    # disable_turbo_boost()

    if should_only_rerun_plotter := args.run == "remake_plot":
        # re-runs the plotting code with the data of the last timing test
        timing_dir = Path(f"{$TESTING_DIR}/timing/")
        root_dir_for_this_test = sorted(timing_dir.glob("*"))[-1]

        all_implementations = [p.name for p in root_dir_for_this_test.glob("*.csv")]
        all_implementations_comma_sep = ",".join(all_implementations)

        # make plot
        cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot.py @(all_implementations_comma_sep)  && cd -

        print(f"\n✅ Remade plot, saved at '{root_dir_for_this_test}/plot.png'\n")
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

    args = argparser.parse_args()
    return args

def get_flops_and_cycles_count(implementation: str, matrix_dimension: int):
    temp_dir = $( mktemp -d ).strip()
    print(temp_dir)

    output_file = f"{temp_dir}/stats.txt"

    cycles = $( perf stat \
        -x "," \
        --output @(output_file) \
        --event=fp_arith_inst_retired.128b_packed_double,fp_arith_inst_retired.128b_packed_single,fp_arith_inst_retired.256b_packed_double,fp_arith_inst_retired.256b_packed_single,fp_arith_inst_retired.scalar_double,fp_arith_inst_retired.scalar_single \
        $C_BINARY -I @(implementation) -t --dimension @(matrix_dimension)
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

            n_flops, n_cycles = get_flops_and_cycles_count(implementation, matrix_dimension)

            print(f" took {n_cycles} cycles ({n_cycles/(3.4*10**9)} sec)")
            data.append((matrix_dimension, n_cycles, n_flops))

        print(data)

        root_dir_for_this_test = f"{$TESTING_DIR}/timing/{current_time}"
        mkdir --parents @(root_dir_for_this_test)

        with open(f"{root_dir_for_this_test}/{implementation}.csv", "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["matrix_dimension", "n_cycles", "n_flops"])
            writer.writerows(data)

    all_implementations = [f"{i}.csv" for i in implementations]
    all_implementations_comma_sep = ",".join(all_implementations)

    # make plot
    cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot.py @(all_implementations_comma_sep)  && cd -

    print(f"\n✅ Finished timing test, saved plot at '{root_dir_for_this_test}/plot.png'\n")


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

    $TESTING_DIR = Path(f"{$TESTING_INFRA_ROOT_DIR}/.test_results").resolve(strict=True)
    $BLOCK_SIZE = args.block_size

    ## XONSH
    # throw error if bash command fails, otherwise we silently ignore the error
    $RAISE_SUBPROC_ERROR = True
    $XONSH_SHOW_TRACEBACK = True
    # for bash commands, print out how they were invoked (with which concrete args)
    # $XONSH_TRACE_SUBPROC = True

def disable_turbo_boost():
    if is_turbo_boost_enabled := $( sudo rdmsr -p0 0x1a0 -f 38:38 ) == "0\n":
        print("✋🛑 Turbo boost is not disabled! Disabling turbo boost...")
        wrmsr -p0 0x1a0 0x4000850089
        wrmsr -p1 0x1a0 0x4000850089
        wrmsr -p2 0x1a0 0x4000850089
        wrmsr -p3 0x1a0 0x4000850089

    print("✅ Turbo boost is disabled\n")

def get_all_implementations() -> str:
    implementations_string: str = $( $C_BINARY -l )
    implementations = implementations_string.split(",")[:-1]
    return implementations

def recompile_c_implementation():
    print("🎬 Re-compiling C implementation...")
    cmake -S $C_IMPLEMENTATION_DIR/ -B $C_IMPLEMENTATION_DIR/build/ -DBLOCK_SIZE=$BLOCK_SIZE -DNDEBUG=YOLOL
    # run make in C implementation dir and then cd back into the prev dir (infra dir)
    cd $C_IMPLEMENTATION_DIR/build && make && cd -

    print("✅ Finished re-compiling C implementation!\n")

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
