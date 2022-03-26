# TomoLearn
A simple C++ library to experiment with the Radon transformation and Filtered Backprojection in parallel beam geometry.

# Dependencies
TomoLearn depends on the following Ubuntu packages:

- `libeigen3-dev`
- `cimg-dev`
- `libfftw3-dev`
- `imagemagick` (needed for png read, could switch to `libpng`)
- `python3-dev`, `python3-numpy`, `python3-matplotlib` (needed for 1D plots; One can use anaconda instead but set corresponding cmake_cache entry)

For building:
- compiler: clang++
- Build system generator: cmake >=3.20 (`sudo snap install cmake --classic`)
- Build system: make  

## CMAKE_CACHE entries
-  `ENABLE_CUDA` bool
-  `ENABLE_SANITIZER_ADDRESS` bool; for line numbers set `ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer` environment variable
-  `ENABLE_SANITIZER_MEMORY` bool
-  `Python3_ROOT_DIR` needed in case of Anaconda Python `/home/szakal/anaconda3/`

## Building and running: 
0. Go to the root directory of the project 
1. `mkdir build`
2. `cd build`
3. `cmake ../ -DCMAKE_BUILD_TYPE=Release`   NOTE: Execution is significantly slower in Debug config 
4. `make`
5. `cd ..`
6. `./build/tomoLearn` 


## Import project to eclipse:
1. Install the cmake4eclipse plugin
2. Follow the instructions in Help -> Help Contents -> CMake4eclipse user guide
3. Set `Unix Makefiles` as the default build system in Windows -> Preferences -> cmake4eclipse
4. Set the compiler to clang by setting CC=clang and CXX=clang++ environment variables in Project Preferences -> C/C++ Build -> Environment

## Select CUDA language settings provider
1. In Preferences -> C/C++ general -> Preprocessor Includes -> Providers -> select "Nvcc Builtins Provider"
2. At the same tab, select "Use global provider shared between projects"
