cmake_minimum_required(VERSION 3.15)

project(TomoLearnLibrary CXX)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)

#Set the available SSE/AVX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# standard compiler warnings
add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# sanitizer options if supported by compiler
include(cmake/Sanitizers.cmake)
enable_sanitizers(project_options)

add_subdirectory(external) #This contains the nontrivial external dependencies
add_subdirectory(tomoLearnLib)
add_subdirectory(tests)
add_subdirectory(benchmarks)

