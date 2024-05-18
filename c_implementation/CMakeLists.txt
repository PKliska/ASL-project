cmake_minimum_required(VERSION 3.9.2)
project(cavity_flow C)

add_subdirectory(external/cargs)


add_executable(cavity_flow src/baseline_simulation.c src/main.c
		           src/simulation.c
			   src/preallocated_simulation.c
			   src/faster_math_simulation.c
			   src/data_trans_simulation.c
			   src/trapeze_simulation.c
			   src/precomputed_trapeze_simulation.c
			   src/utils.c)

target_compile_options(cavity_flow PUBLIC -g -O3 -Wall -Wextra -mavx2 -mfma)
target_include_directories(cavity_flow PUBLIC include)
target_link_libraries(cavity_flow cargs)
set_target_properties(cavity_flow
		      PROPERTIES
		      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
