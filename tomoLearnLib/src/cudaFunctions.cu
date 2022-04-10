#include <Gen1CT.hpp>

#include <config.h>
#include <cstdio>

#if ENABLE_CUDA

#include <cuda_runtime.h>
//#include <device_launch_parameters.h>

#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

template<typename T>
void check(T err, const char* const func, const char* const file, const int line) {
  if (err != cudaSuccess) {
    std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
    std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
    exit(1);
  }
}


void Gen1CT::printGpuParameters(){
	int nDevices;

	cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++) {
	    cudaDeviceProp prop;
	    cudaGetDeviceProperties(&prop, i);
	    printf("Device Number: %d\n", i);
	    printf("  Device name: %s\n", prop.name);
	    printf("  Memory Clock Rate (KHz): %d\n",
	           prop.memoryClockRate);
	    printf("  Memory Bus Width (bits): %d\n",
	           prop.memoryBusWidth);
	    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
	           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
	}
}

/***
 * Ray driven projection (Hao Gao method) implemented on the GPU
 * @param actualPhantom The phantom which should be projected
 * @param angles Vector of projection angles
 * @return Returns the sinogram of the actualPhantom
 */
Eigen::MatrixXd Gen1CT::project_rayDriven_GPU(const Phantom& actualPhantom,
		                    const Eigen::VectorXd& angles){
	std::cout << std::endl << "Projection with ray-driven method started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(static_cast<long>(pixNum), static_cast<long>(numAngles));

	auto pixSizes = actualPhantom.getPixSizes();
	auto numberOfPixels = actualPhantom.getNumberOfPixels();

	////////////////////////////////////////////////////////////////////////////////
	/// Projection with ray-driven method on GPU STARTS here !!
	////////////////////////////////////////////////////////////////////////////////
	double *d_phantom;
	std::array<int, 2> numPixels = actualPhantom.getNumberOfPixels();
	checkCudaErrors(cudaMalloc(&d_phantom, sizeof(double) * numPixels[0]*numPixels[1] ));
	checkCudaErrors(cudaMemcpy(d_phantom, actualPhantom.getDataAsEigenMatrixRef().data(),
			                   sizeof(double) * numPixels[0]*numPixels[1], cudaMemcpyHostToDevice));



	////////////////////////////////////////////////////////////////////////////////
	/// Projection with ray-driven method on GPU ENDS here !!
	////////////////////////////////////////////////////////////////////////////////

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
	std::cout << "Projection with ray-driven method took " << duration.count() << " milliseconds" << std::endl;

	return sinogram;
}

#endif
