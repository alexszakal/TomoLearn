cmake_minimum_required(VERSION 3.15)

project(TomoLearn CXX)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)

# standard compiler warnings
add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# sanitizer options if supported by compiler
include(cmake/Sanitizers.cmake)
enable_sanitizers(project_options)

add_subdirectory(TomoLearnS)

add_executable(tomoLearn "main.cpp")
target_link_libraries(tomoLearn PUBLIC TomoLearnS project_options)
