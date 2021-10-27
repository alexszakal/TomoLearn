add_library(tomoLearnLib STATIC src/CTScan.cpp 
                                src/Filter.cpp 
                                src/Gen1CT.cpp 
                                src/Object2D.cpp 
                                src/Phantom.cpp)

option(ENABLE_CUDA "Turn on CUDA reconstructions" ON)
  
#configure_file(config.h.in config.h)
#target_include_directories(tomoLearnLib PUBLIC "${TomoLearnS_BINARY_DIR}")

target_compile_features(tomoLearnLib PRIVATE cxx_std_17)

configure_file(src/config.h.in config.h)
message(STATUS "config.h written to: "  "${CMAKE_CURRENT_BINARY_DIR}")
target_include_directories(tomoLearnLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(tomoLearnLib PUBLIC include)
target_include_directories(tomoLearnLib PRIVATE src)

target_include_directories(tomoLearnLib SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../matplotlibcpp)
target_include_directories(tomoLearnLib SYSTEM PUBLIC /usr/include/eigen3)
target_include_directories(tomoLearnLib SYSTEM PUBLIC /usr/include/python2.7)

target_link_libraries(tomoLearnLib PUBLIC python2.7)
target_link_libraries(tomoLearnLib PUBLIC pthread)
target_link_libraries(tomoLearnLib PUBLIC fftw3)
target_link_libraries(tomoLearnLib PUBLIC m)
target_link_libraries(tomoLearnLib PUBLIC X11)

target_link_libraries(tomoLearnLib PUBLIC project_options project_warnings)



