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
./run-script.sh ./testing_infra.xsh
```

