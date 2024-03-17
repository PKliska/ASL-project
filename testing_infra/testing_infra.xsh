import csv
from typing import Any

def main():
    # throw error if bash command fails, otherwise we silently ignore the error
    $RAISE_SUBPROC_ERROR = True
    # for bash commands, print out how they were invoked (with which concrete args)
    # $XONSH_TRACE_SUBPROC = True

    print("\n\n")

    $TESTING_DIR = "./.test_results"
    mkdir --parents "$TESTING_DIR"

    make ../c_implementation > /dev/null # compile C code

    run_consystency_test()

def run_consystency_test():
    for i in range(1000):
        for _ in range(2):
            print(f"n_simulation_iterations = {i} ", end="")
            check_if_c_outputs_consistent_for(i)

def check_if_c_outputs_consistent_for(n_simulation_iterations: int = 100):
    """Checks if the C implementation produces the same output CSV file if given the
    same arguments."""

    c_output_path_1 = f"{$TESTING_DIR}/output_c_1.csv"
    c_output_path_2 = f"{$TESTING_DIR}/output_c_2.csv"


    ../c_implementation/bin/baseline --output_file=@(c_output_path_1) --num_iter=@(n_simulation_iterations)
    ../c_implementation/bin/baseline --output_file=@(c_output_path_2) --num_iter=@(n_simulation_iterations)

    compare_files_command = !( cmp --silent -- @(c_output_path_1) @(c_output_path_2 ))
    if are_files_different := compare_files_command.returncode != 0:
        print("❌")
        raise Exception("Running the C implementation with the same inputs produced different outputs!")
    else:
        print("✅")

def check_if_c_output_matches_python():
    n_simulation_iterations = 100

    # TODO: run python implementation & save output


    # compile and run C implementation & save output
    c_output_path = f"{$TESTING_DIR}/output_c.csv"

    make
    ../c_implementation/bin/baseline --output_file=@(c_output_path) --num_iter=@(n_simulation_iterations)

    # print(read_matrices_from_csv(c_output_path)[0])

    # TODO: compare the outputs


## HELPER FUNCTIONS
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

## EXECUTE MAIN
main()