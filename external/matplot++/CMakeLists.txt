cmake_minimum_required( VERSION 3.15 )

find_package(Matplot++ QUIET)

if(Matplot++_FOUND)
	message(STATUS "Matplot++ was found system-wide")
else()
	message(STATUS "Matplot++ was not found system-wide, downloading from github")
	
	include(FetchContent)
	
	FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master)
        
    FetchContent_GetProperties(matplotplusplus)
	if(NOT matplotplusplus_POPULATED)
    	FetchContent_Populate(matplotplusplus)
    	message(STATUS "Matplot++ sourceDir: " ${matplotplusplus_SOURCE_DIR})
    	set(BUILD_SHARED_LIBS ON) #matplotplusplus builds static libs by default
    	add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
    
endif()