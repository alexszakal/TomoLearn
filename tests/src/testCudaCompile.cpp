#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>

#include <CImg.h>

#ifdef Success       //Because otherwise Eigen not compile (EIgen <-> CImg interference)
  #undef Success
#endif

#include <Object2D.hpp>
#include <Gen1CT.hpp>

#include <config.h>

void testCudaCompile();

int main(){

#if ENABLE_CUDA
	std::cout << "\n \n CUDA enabled!!!!" ;
	testCudaCompile();
#else
	std::cout << "\n \n CUDA disabled!!!" ;
#endif

	std::cout << "\nPress ENTER to continue!";
	std::cin.get();

	return 0;
}

#if ENABLE_CUDA
void testCudaCompile(){

	std::cout << "\nTest if Cuda device query works" << std::endl;

	int detWidthInMM { 110 };
	int detPixNum { 512 };
	Gen1CT ct(detWidthInMM, detPixNum);

	ct.printGpuParameters();
}
#endif


