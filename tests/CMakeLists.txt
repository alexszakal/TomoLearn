
add_executable(testFBP src/testFBP.cpp)
target_link_libraries(testFBP PRIVATE tomoLearnLib)

add_executable(testMLEM src/testMLEM.cpp)
target_link_libraries(testMLEM PRIVATE tomoLearnLib)

add_executable(testRadonTransform src/testRadonTransform.cpp)
target_link_libraries(testRadonTransform PRIVATE tomoLearnLib)

add_executable(testCudaCompile src/testCudaCompile.cpp)
target_link_libraries(testCudaCompile PRIVATE tomoLearnLib)

add_executable(testRayDrivenProj src/testRayDrivenProj.cpp)
target_link_libraries(testRayDrivenProj PRIVATE tomoLearnLib)

add_executable(testSPS src/testSPS.cpp)
target_link_libraries(testSPS PRIVATE tomoLearnLib)

add_executable(testReconstructions src/testReconstructions.cpp)
target_link_libraries(testReconstructions PRIVATE tomoLearnLib)