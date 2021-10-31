
add_executable(testFBP src/testFBP.cpp)
target_link_libraries(testFBP PRIVATE tomoLearnLib)

add_executable(testRadonTransform src/testRadonTransform.cpp)
target_link_libraries(testRadonTransform PRIVATE tomoLearnLib)

add_executable(testCudaCompile src/testCudaCompile.cpp)
target_link_libraries(testCudaCompile PRIVATE tomoLearnLib)

