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

# runs just the (short) timing test to generate the plot, quick but only times the algo for a few smaller sizes
./run-script.sh ./testing_infra.xsh --run timing --short-test --implementation baseline

# runs the long timing test which takes way longer to complete, but tests more sizes (like >10)
./run-script.sh ./testing_infra.xsh --run timing --long-test

# runs the long timing test for multiple implementations
./run-script.sh ./testing_infra.xsh --run timing --long-test --implementation baseline preallocated


```

