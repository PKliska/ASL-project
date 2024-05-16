import csv
from typing import Any
import sys
from time import strftime
from pathlib import Path
from argparse import ArgumentParser


def main():
    disable_turbo_boost()

    args = parse_cli_args()
    set_global_variables(args)

    recompile_c_implementation()
    mkdir --parents "$TESTING_DIR"

    # based on arguments, run timing plot, consystency, correctness, etc.
    if args.run == "timing":
        if args.long_test:
            run_timing_test(100, 1901, 200)
        else:
            run_timing_test(100, 301, 100)
    elif args.run == "correctness":
        run_correctness_test()
    elif args.run == "consystency":
        run_consystency_test()
    elif args.run == "all":
        run_correctness_test()
        run_consystency_test()
        if args.long_test:
            run_timing_test(100, 1901, 200)
        else:
            run_timing_test(100, 301, 100)

    else:
        raise Exception(f"Can't run action '{args.run}', invalid option. Did you mispell?")

def parse_cli_args():
    argparser = ArgumentParser(description="")
    argparser.add_argument("--run", metavar="ACTION", default="all", help="Action to execute (e.g. run correctness tests)")
    argparser.add_argument("--implementation", default="baseline", help="Chose which implementation to run")
    group = argparser.add_mutually_exclusive_group()
    group.add_argument("--short-test", action="store_true", default=True, help="Quickly gets some (rough) results")
    group.add_argument("--long-test", action="store_true", help="Slowly gets detailed results, runs more iterations of the algo")

    args = argparser.parse_args()
    return args

def run_timing_test(start: int, stop: int, step: int):
    print("\n🟠 Starting timing test...")

    current_time = strftime("%Y-%m-%d %H:%M:%S")

    for implementation in get_all_implementations():
        data = []
        for matrix_dimension in range(start, stop, step):
            output = $( $C_BINARY -I @(implementation) -t --dimension @(matrix_dimension) )
            n_cycles = float(output.strip().split()[-1])

            print(f"Simulating {implementation} with dimension={matrix_dimension} took {n_cycles} cycles ({n_cycles/(2.5*10**9)} sec)")
            data.append((matrix_dimension, n_cycles))

        print(data)

        root_dir_for_this_test = f"{$TESTING_DIR}/timing/{current_time}"
        mkdir --parents @(root_dir_for_this_test)

        with open(f"{root_dir_for_this_test}/{implementation}.csv", "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["matrix_dimension", "n_cycles"])
            writer.writerows(data)

    all_implementations = [f"{i}.csv" for i in get_all_implementations()]
    all_implementations_comma_sep = ",".join(all_implementations)

    # make plot
    cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot.py @(all_implementations_comma_sep)  && cd -

    print(f"\n✅ Finished timing test, saved plot at '{root_dir_for_this_test}/plot.png'\n")


def run_correctness_test():
    print("\n🟠 Starting correctness test...")

    # test for 4, 5 (for 1..3 it fails cuz of a formatting mismatch...), some regularly spaced ones and big ones just in case
    dimensions_which_to_test = [4, 5, 289, 800, 950] + list(range(10, 123, 7))
    dimensions_which_to_test = sorted(list(set(dimensions_which_to_test))) # remove duplicates & sort

    is_some_result_incorrect = False
    for i in dimensions_which_to_test:
        print(f"matrix_dimension = {i} ", end="", flush=True)
        try:
            check_if_c_output_matches_python_output_for(i)
        except Exception:
            is_some_result_incorrect = True
            pass

    # manual test for nan/inf
    very_large_matrix_dimension = 1900
    print(f"matrix_dimension (just check for NaN/Inf, don't compare to python impl.) = {very_large_matrix_dimension} ", end="", flush=True)

    root_dir_for_very_large_matrix = f"{$TESTING_DIR}/correctness/n_{very_large_matrix_dimension}"
    mkdir --parents @(root_dir_for_very_large_matrix)

    c_output_path = f"{root_dir_for_very_large_matrix}/output_c_{$IMPLEMENTATION}.csv"
    $C_BINARY  -I "$IMPLEMENTATION" --output_file=@(c_output_path) --dimension=@(very_large_matrix_dimension)

    command = !( grep -q -Ei "nan|inf" @(c_output_path) )

    if does_matrix_contain_nan_or_inf := not has_command_failed(command):
        print("❌      ❗️ contains NaN or Inf")
        raise Exception("❗️❗️ Correctness tests failed!")
    else:
        print("✅")


    if is_some_result_incorrect:
        raise Exception("❗️❗️ Correctness tests failed!")
    print(f"\n✅ Finished correctness test!\n")


def run_consystency_test():
    print(f"\n🟠 Starting consystency test for '{$IMPLEMENTATION}'...")

    for i in range(10, 101, 10):
        print(f"matrix_dimension = {i} ", end="", flush=True)
        check_if_c_outputs_consistent_for(i)

    print(f"\n✅ Finished consystency test!\n")


def check_if_c_outputs_consistent_for(matrix_dimension):
    """Checks if the C implementation produces the same output CSV file if given the
    same arguments."""

    mkdir -p "$TESTING_DIR/consystency/"
    c_output_path_1 = f"{$TESTING_DIR}/consystency/output_c_1_{$IMPLEMENTATION}.csv"
    c_output_path_2 = f"{$TESTING_DIR}/consystency/output_c_2_{$IMPLEMENTATION}.csv"

    $C_BINARY -I "$IMPLEMENTATION" --output_file=@(c_output_path_1) --dimension=@(matrix_dimension)
    $C_BINARY -I "$IMPLEMENTATION" --output_file=@(c_output_path_2) --dimension=@(matrix_dimension)

    compare_files_command = !( cmp --silent -- @(c_output_path_1) @(c_output_path_2 ))
    if are_files_different := has_command_failed(compare_files_command):
        print("❌")
        raise Exception("Running the C implementation with the same inputs produced different outputs!")
    else:
        print("✅")

def check_if_c_output_matches_python_output_for(matrix_dimension):
    root_dir_for_this_test = f"{$TESTING_DIR}/correctness/n_{matrix_dimension}"
    mkdir --parents @(root_dir_for_this_test)

    # run python implementation & save output
    python_output_path = f"{root_dir_for_this_test}/output_python.csv"
    python ./python_implementation.py @(python_output_path) @(matrix_dimension)

    # run C implementation & save output
    c_output_path = f"{root_dir_for_this_test}/output_c_{$IMPLEMENTATION}.csv"
    $C_BINARY  -I "$IMPLEMENTATION" --output_file=@(c_output_path) --dimension=@(matrix_dimension)

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
    $IMPLEMENTATION = args.implementation

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
    cmake -S $C_IMPLEMENTATION_DIR/ -B $C_IMPLEMENTATION_DIR/build/
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
