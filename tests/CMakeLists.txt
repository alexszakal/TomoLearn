
add_executable(testFBP src/testFBP.cpp)
target_link_libraries(testFBP PRIVATE tomoLearnLib)

add_executable(testMLEM src/testMLEM.cpp)
target_link_libraries(testMLEM PRIVATE tomoLearnLib)

add_executable(testRadonTransform src/testRadonTransform.cpp)
target_link_libraries(testRadonTransform PRIVATE tomoLearnLib)

add_executable(testCudaCompile src/testCudaCompile.cpp)
target_link_libraries(testCudaCompile PRIVATE tomoLearnLib)

add_executable(testHaoGaoCPU src/testHaoGaoCPU.cpp)
target_link_libraries(testHaoGaoCPU PRIVATE tomoLearnLib)

