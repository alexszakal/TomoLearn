cmake_minimum_required(VERSION 3.20)

message(STATUS "Searching for Eigen3  library...")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set_target_properties(Eigen3::Eigen PROPERTIES IMPORTED_GLOBAL TRUE)
