import csv
from typing import Any
import sys
from time import strftime
from pathlib import Path
from argparse import ArgumentParser


def main():
    args = parse_cli_args()
    set_global_variables()

    recompile_c_implementation()
    mkdir --parents "$TESTING_DIR"

    # based on arguments, run timing plot, consystency, correctness, etc.
    if args.run == "timing":
        if args.quick:
            run_timing_test(100, 5000, 2000)
        else:
            run_timing_test(100, 50000, 10000)
    elif args.run == "correctness":
        run_correctness_test()
    elif args.run == "consystency":
        run_consystency_test()
    elif args.run == "all":
        run_correctness_test()
        run_consystency_test()
        run_timing_test(100, 50000, 10000)
    else:
        raise Exception(f"Can't run action '{args.run}', invalid option. Did you mispell?")

def parse_cli_args():
    argparser = ArgumentParser(description="")
    argparser.add_argument("--run", metavar="ACTION", default="all", help="Action to execute (e.g. run correctness tests)")
    group = argparser.add_mutually_exclusive_group()
    group.add_argument("--quick", action="store_true", default=True, help="Quickly gets some (rough) results")
    group.add_argument("--slow", action="store_true", help="Slowly gets detailed results, runs more iterations of the algo")

    args = argparser.parse_args()
    return args

def run_timing_test(start: int, stop: int, step: int):
    print("\n🟠 Starting timing test...")

    current_time = strftime("%Y-%m-%d %H:%M:%S")

    for implementation in ["baseline", "preallocated"]:
        data = []
        for n_simulation_iterations in range(start, stop, step):
            output = $( $C_IMPLEMENTATION_DIR/build/bin/cavity_flow -I @(implementation) -t --num_iter @(n_simulation_iterations) )
            n_cycles = float(output.strip().split()[-1])

            print(f"Simulating {implementation} with n_iter={n_simulation_iterations} took {n_cycles} cycles ({n_cycles/(2.5*10**9)} sec)")
            data.append((n_simulation_iterations, n_cycles))

        print(data)

        root_dir_for_this_test = f"{$TESTING_DIR}/timing/{current_time}"
        mkdir --parents @(root_dir_for_this_test)

        with open(f"{root_dir_for_this_test}/{implementation}.csv", "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["n_simulation_iterations", "n_cycles"])
            writer.writerows(data)

    # make plot
    cd @(root_dir_for_this_test) && python $TESTING_INFRA_ROOT_DIR/timing_plot.py && cd -

    print(f"\n✅ Finished timing test, saved plot at '{root_dir_for_this_test}/plot.png'\n")


def run_correctness_test():
    print("\n🟠 Starting correctness test...")

    iterations_which_to_test = [0, 1] + list(range(0, 1000, 77)) # test for 0, 1 & some odd ones just in case
    iterations_which_to_test = sorted(list(set(iterations_which_to_test))) # remove duplicates & sort

    is_some_result_incorrect = False
    for i in iterations_which_to_test:
        print(f"n_simulation_iterations = {i} ", end="")
        try:
            check_if_c_output_matches_python_output_for(i)
        except Exception:
            is_some_result_incorrect = True
            pass

    if is_some_result_incorrect:
        raise Exception("❗️❗️ Correctness tests failes!")
    print(f"\n✅ Finished correctness test!\n")


def run_consystency_test():
    print("\n🟠 Starting consystency test...")

    for i in range(0, 1000, 77):
        print(f"n_simulation_iterations = {i} ", end="")
        check_if_c_outputs_consistent_for(i)

    print(f"\n✅ Finished consystency test!\n")


def check_if_c_outputs_consistent_for(n_simulation_iterations: int = 100):
    """Checks if the C implementation produces the same output CSV file if given the
    same arguments."""

    c_output_path_1 = f"{$TESTING_DIR}/output_c_1.csv"
    c_output_path_2 = f"{$TESTING_DIR}/output_c_2.csv"

    $C_IMPLEMENTATION_DIR/build/bin/cavity_flow -I "baseline" --output_file=@(c_output_path_1) --num_iter=@(n_simulation_iterations)
    $C_IMPLEMENTATION_DIR/build/bin/cavity_flow -I "baseline" --output_file=@(c_output_path_2) --num_iter=@(n_simulation_iterations)

    compare_files_command = !( cmp --silent -- @(c_output_path_1) @(c_output_path_2 ))
    if are_files_different := has_command_failed(compare_files_command):
        print("❌")
        raise Exception("Running the C implementation with the same inputs produced different outputs!")
    else:
        print("✅")

def check_if_c_output_matches_python_output_for(n_simulation_iterations: int = 100):
    root_dir_for_this_test = f"{$TESTING_DIR}/correctness/n_{n_simulation_iterations}"
    mkdir --parents @(root_dir_for_this_test)

    # run python implementation & save output
    python_output_path = f"{root_dir_for_this_test}/output_python.csv"
    python ./python_implementation.py @(python_output_path) @(n_simulation_iterations)

    # run C implementation & save output
    c_output_path = f"{root_dir_for_this_test}/output_c.csv"
    $C_IMPLEMENTATION_DIR/build/bin/cavity_flow --output_file=@(c_output_path) --num_iter=@(n_simulation_iterations)

    # compare the outputs
    compare_files_command = !( cmp --silent -- @(c_output_path) @(python_output_path) )
    if are_files_different := has_command_failed(compare_files_command):
        print("❌")
        raise Exception("Error: C implementation produced different result compared to Python implementation!")
    else:
        print("✅")


## HELPER FUNCTIONS
def set_global_variables():
    ## VARS FOR THIS PROGRAM
    $TESTING_INFRA_ROOT_DIR = Path(__file__).parent.resolve(strict=True)
    $C_IMPLEMENTATION_DIR = Path(f"{$TESTING_INFRA_ROOT_DIR}/../c_implementation").resolve(strict=True)

    $TESTING_DIR = Path(f"{$TESTING_INFRA_ROOT_DIR}/.test_results").resolve(strict=True)

    ## XONSH
    # throw error if bash command fails, otherwise we silently ignore the error
    $RAISE_SUBPROC_ERROR = True
    $XONSH_SHOW_TRACEBACK = True
    # for bash commands, print out how they were invoked (with which concrete args)
    # $XONSH_TRACE_SUBPROC = True

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
