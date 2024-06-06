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
			   src/skewed_simulation.c
			   src/uskewed_simulation.c
			   src/subskewed_simulation.c
			   src/utils.c)

add_custom_command(
	OUTPUT include/generated_simulation.h src/generated_simulation.c
	COMMAND python ${PROJECT_SOURCE_DIR}/codegen/codegen.py ${PROJECT_BINARY_DIR} 4 2
	DEPENDS codegen/codegen.py
)
target_include_directories(cavity_flow PRIVATE ${PROJECT_BINARY_DIR}/include)
target_sources(cavity_flow PRIVATE ${PROJECT_BINARY_DIR}/src/generated_simulation.c)



target_compile_options(cavity_flow PUBLIC -g -fopt-info -O3 -Wall -Wextra -mavx2 -mfma -ffast-math)
target_include_directories(cavity_flow PUBLIC include)
target_link_libraries(cavity_flow cargs)

# add new preprocessor options here
set(preprocessor_options
	NDEBUG

	TRAPEZE_CUTOFF_VOLUME
	TRAPEZE_CUTOFF_AREA

	SKEWING_BLOCK_SIZE_X
	SKEWING_BLOCK_SIZE_Y
	SKEWING_TIMESTEPS

	SSKEWING_BLOCK_SIZE_X
	SSKEWING_BLOCK_SIZE_Y
	SSKEWING_SUBBLOCK_SIZE_X
	SSKEWING_SUBBLOCK_SIZE_Y
	SSKEWING_TIMESTEPS

	USKEWING_BLOCK_SIZE_X
	USKEWING_BLOCK_SIZE_Y
	USKEWING_TIMESTEPS
)

foreach(option IN LISTS preprocessor_options)
	if(DEFINED ${option})
		message(STATUS ${option}=${${option}})
		target_compile_definitions(cavity_flow PUBLIC
			${option}=${${option}})
	endif()
endforeach()
set_target_properties(cavity_flow
		      PROPERTIES
		      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
