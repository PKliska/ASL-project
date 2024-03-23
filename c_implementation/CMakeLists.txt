cmake_minimum_required(VERSION 3.9.2)
project(cavity_flow C)

add_subdirectory(external/cargs)

add_executable(cavity_flow src/baseline_simulation.c src/main.c)

target_compile_options(cavity_flow PUBLIC -O3 -Wall -Wextra)
target_include_directories(cavity_flow PUBLIC include)
target_link_libraries(cavity_flow cargs)
set_target_properties(cavity_flow
		      PROPERTIES
		      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
