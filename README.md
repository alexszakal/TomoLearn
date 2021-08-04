# TomoLearn
A simple C++ library to experiment with the Radon transformation and Filtered Backprojection in parallel beam geometry.

# Dependencies
TomoLearn depends on the following Ubuntu packages:
- CImg
- fftw3
- python2.7, python-numpy, python-matplotlib  (for 1D plotting)
- compiler: clang++ 

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
3. Set the compiler to clang by setting CC=clang and CXX=clang++ environment variables in Project Preferences -> C/C++ Build -> Environment

## Select CUDA language settings provider
1. In Preferences -> C/C++ general -> Preprocessor Includes -> Providers -> select "Nvcc Builtins Provider"
2. At the same tab, select "Use global provider shared between projects"
