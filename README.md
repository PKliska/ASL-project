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
./run-script.sh ./testing_infra.xsh # runs all tests (correctness, consystency & timing) for baseline

# choose which implementation to test, if not implementation passed "baseline" is used
./run-script.sh ./testing_infra.xsh --implementation preallocated

# runs just the (short) timing test to generate the plot, quick but only times the algo for 3 "sizes"
./run-script.sh ./testing_infra.xsh --run timing --short

# runs the long timing test which takes way longer to complete, but tests more sizes (like >10)
./run-script.sh ./testing_infra.xsh --run timing --long

# runs correctness test
./run-script.sh ./testing_infra.xsh --run correctness
```

