# Cavity flow simulation (team 22)

An Advanced Systems Lab project written by: Selin Baraş, Peter Horcic,
Filip Jakšić, and Pavel Kliska.

## Building and running the C implementation
```bash
cd c_implementation
mkdir build
cd build
cmake ..
make
bin/cavity_flow
```
To do a debug build add the `-DCMAKE_BUILD_TYPE=Debug` flag to the cmake
command.

## Running the testing infrastructure
```bash
cd testing_infra

# create virtual env as not to pollute your global python install
python3 -m venv .venv
source .venv/bin/activate

# install dependencies
python3 -m pip install -r requirements.txt

# actually run the testing infrastructure
./run-script.sh ./testing_infra.xsh # runs all tests (correctness, consystency & timing)

# run the testing infrastructure for a specific implementation, runs all tests (correctness, consystency & timing)
./run-script.sh ./testing_infra.xsh --implementation baseline

# run all tests for multiple implementations (e.g. for baseline & prealloc)
./run-script.sh ./testing_infra.xsh --implementation baseline preallocated

# run all tests for multiple implementations (e.g. for baseline & prealloc) for matrices of dimensions [10, 20, 30]
./run-script.sh ./testing_infra.xsh --implementation baseline preallocated --matrix-dimensions 10 20 30

# runs just the timing test to generate the plot
./run-script.sh ./testing_infra.xsh --run timing --implementation baseline

# runs the timing test for multiple implementations
./run-script.sh ./testing_infra.xsh --run timing --implementation baseline preallocated

# runs the timing test for multiple implementations and for specific matrix dimensions
./run-script.sh ./testing_infra.xsh --run timing --implementation baseline preallocated --matrix-dimensions 32 64 96

# Remake the plot for the last timing test, i.e. doesn't run anything, just takes the measurements in that folder and
# re-generates the plot. Can be useful if you wanna combine multiple previous timing tests into 1 plot, e.g. you run a
# separate timing test just to baseline and another separata for prealloc, then you can just copy/paste the CSV files
# from one test into the other one and generate a plot which contains both, instead of wasting time by running a
# separate timing test which contains both :)
./run-script.sh ./testing_infra.xsh --run remake_plot


```

