cmake_minimum_required( VERSION 3.15)

add_library(matplotlib_cpp INTERFACE)
set_property(TARGET matplotlib_cpp PROPERTY CXX_STANDARD 11)

message(STATUS "Adding the matplotlibcpp header-only library")

find_package(Python3 COMPONENTS Interpreter Development NumPy )
target_include_directories(matplotlib_cpp INTERFACE ${Python3_INCLUDE_DIRS})
target_link_libraries(matplotlib_cpp INTERFACE Python3::Python Python3::Module Python3::NumPy)

message(STATUS "python3 interpreter ID: " ${Python3_INTERPRETER_ID})
message(STATUS "python3 include dirs: " ${Python3_INCLUDE_DIRS})

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
    matplotlibcpp
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        ef0383f1315d32e0156335e10b82e90b334f6d9f 
)
FetchContent_Populate(matplotlibcpp)
message(STATUS "matplotlibSourceDir: " ${matplotlibcpp_SOURCE_DIR})
target_include_directories(matplotlib_cpp SYSTEM INTERFACE ${matplotlibcpp_SOURCE_DIR})
