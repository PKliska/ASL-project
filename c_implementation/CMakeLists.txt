cmake_minimum_required(VERSION 3.9.2)
project(cavity_flow C)

add_subdirectory(external/cargs)

# Create the main executable target
add_executable(cavity_flow src/baseline_simulation.c src/main.c
			   src/simulation.c
			   src/preallocated_simulation.c
			   src/faster_math_simulation.c
			   src/fasterfaster_math_simulation.c
			   src/data_trans_simulation.c
			   src/vanilla_avx_simulation.c
			   src/aligned_avx_simulation.c
			   src/trapeze_simulation.c
			   src/precomputed_trapeze_simulation.c
			   src/cutoff_trapeze_simulation.c
			   src/blocking_simulation.c
			   src/utils.c)

# Define a cache variable for the custom flag and define value
set(BLOCK_SIZE "" CACHE STRING "Custom define value")

# Set compilation options and include directories for the main target
if (BLOCK_SIZE)
	message(STATUS "❗️ Using BLOCK SIZE: ${BLOCK_SIZE}")
	add_compile_definitions(BLOCK_SIZE=${BLOCK_SIZE})
	target_compile_options(cavity_flow PUBLIC -g -O3 -Wall -Wextra -mavx2 -mfma)
else()
	message(STATUS "❗️ Using BLOCK SIZE: 8")
	add_compile_definitions(BLOCK_SIZE=8)
	target_compile_options(cavity_flow PUBLIC -g -O3 -Wall -Wextra -mavx2 -mfma)
endif()

target_include_directories(cavity_flow PUBLIC include)
target_link_libraries(cavity_flow cargs)
set_target_properties(cavity_flow
		      PROPERTIES
		      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
