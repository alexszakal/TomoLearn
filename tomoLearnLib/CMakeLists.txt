cmake_minimum_required(VERSION 3.20)

add_library(tomoLearnLib SHARED src/CTScan.cpp 
                                src/Filter.cpp 
                                src/Gen1CT.cpp 
                                src/Object2D.cpp 
                                src/Phantom.cpp)

option(ENABLE_CUDA "Turn on CUDA reconstructions" ON)

if(ENABLE_CUDA)
    enable_language(CUDA)
    
    find_package(CUDAToolkit REQUIRED)
    if(NOT CUDAToolkitFOUND)
	message(STATUS "CUDA TOOLKIT FOUND")
    else()
	message(WARNING "CUDA TOOLKIT WAS NOT FOUND)")
    endif()

    target_sources(tomoLearnLib PRIVATE src/cudaFunctions.cu)

    set_property(TARGET tomoLearnLib PROPERTY CUDA_ARCHITECTURES OFF)
    message(STATUS "Cmake include directories: " ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_include_directories(tomoLearnLib SYSTEM PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(tomoLearnLib PUBLIC CUDA::cudart)
    set_target_properties(tomoLearnLib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

target_compile_features(tomoLearnLib PRIVATE cxx_std_17)

#Write CMAKE variables to C++ header file
configure_file(src/config.h.in config.h)
message(STATUS "config.h written to: "  "${CMAKE_CURRENT_BINARY_DIR}")
target_include_directories(tomoLearnLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(tomoLearnLib PUBLIC include)
target_include_directories(tomoLearnLib PRIVATE src)

target_link_libraries(tomoLearnLib PUBLIC Eigen3::Eigen)
target_link_libraries(tomoLearnLib PUBLIC matplotlib_cpp)

target_link_libraries(tomoLearnLib PUBLIC pthread)
target_link_libraries(tomoLearnLib PUBLIC fftw3)
target_link_libraries(tomoLearnLib PUBLIC m)
target_link_libraries(tomoLearnLib PUBLIC X11)

target_link_libraries(tomoLearnLib PUBLIC project_options project_warnings)



