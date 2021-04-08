cmake_minimum_required(VERSION 3.10)

project(TomoLearn)

add_subdirectory(TomoLearnS)

add_executable(tomoLearn "main.cpp")
target_compile_features(tomoLearn PRIVATE cxx_std_17)
target_link_libraries(tomoLearn PUBLIC TomoLearnS)


