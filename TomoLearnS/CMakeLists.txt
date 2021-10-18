add_library(TomoLearnS STATIC CTScan.cpp Filter.cpp Gen1CT.cpp Object2D.cpp Phantom.cpp)

target_compile_features(TomoLearnS PRIVATE cxx_std_17)

target_include_directories(TomoLearnS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(TomoLearnS SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../matplotlibcpp)
target_include_directories(TomoLearnS SYSTEM PUBLIC /usr/include/eigen3)
target_include_directories(TomoLearnS SYSTEM PUBLIC /usr/include/python2.7)

target_link_libraries(TomoLearnS PUBLIC python2.7)
target_link_libraries(TomoLearnS PUBLIC pthread)
target_link_libraries(TomoLearnS PUBLIC fftw3)
target_link_libraries(TomoLearnS PUBLIC m)
target_link_libraries(TomoLearnS PUBLIC X11)

target_link_libraries(TomoLearnS PUBLIC project_options project_warnings)



