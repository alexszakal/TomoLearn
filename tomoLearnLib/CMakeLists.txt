cmake_minimum_required(VERSION 3.20)

add_library(tomoLearnLib SHARED src/CTScan.cpp 
                                src/Filter.cpp 
                                src/Gen1CT.cpp 
                                src/Object2D.cpp 
                                src/Phantom.cpp
                                src/Reconst.cpp)
              
option(ENABLE_CUDA "Turn on CUDA reconstructions" OFF)

if(ENABLE_CUDA)
    enable_language(CUDA)
    
    find_package(CUDAToolkit REQUIRED)
    if(NOT CUDAToolkitFOUND)
	message(STATUS "CUDA TOOLKIT FOUND")
    else()
	message(WARNING "CUDA TOOLKIT WAS NOT FOUND)")
    endif()
    
    #target_compile_options(tomoLearnLib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:"">) #Specific options for nvcc 
    
    IF(CMAKE_BUILD_TYPE)  #Check if CMAKE_BUILD_TYPE is defined
    	IF(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    		set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")   #Needed for the CUDA memcheck tool to show line numbers
    	ENDIF()
    ENDIF()
	
    target_sources(tomoLearnLib PRIVATE src/cudaFunctions.cu)

    set_property(TARGET tomoLearnLib PROPERTY CUDA_ARCHITECTURES all-major)
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
target_link_libraries(tomoLearnLib PUBLIC matplot)

target_link_libraries(tomoLearnLib PUBLIC pthread)
target_link_libraries(tomoLearnLib PUBLIC fftw3)
target_link_libraries(tomoLearnLib PUBLIC m)
target_link_libraries(tomoLearnLib PUBLIC X11)

target_link_libraries(tomoLearnLib PUBLIC project_options project_warnings)



