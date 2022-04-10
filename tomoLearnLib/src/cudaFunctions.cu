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

/***
 * Ray driven projection (Hao Gao method) KERNEL
 * @param phantom Pointer to the phantom data
 * @param sinogram Pointer to the resulting sinogram
 * @return Returns void
 */
__global__
void rayDrivenKernel(double* phantom, int numberOfPixelsX, int numberOfPixelsY, double pixSizesX, double pixSizesY,
		             double* sinogram, int numAngles, double* angles, int numDetPixels, double detPixSize ){

	int detPixIdx = blockIdx.x*blockDim.x + threadIdx.x;
	int angIdx =    blockIdx.y*blockDim.y + threadIdx.y;

	//printf("\n detPixIdx: %d, angIdx: %d", detPixIdx, angIdx);

	if ( detPixIdx >= numDetPixels){
		return;
	}
	if (angIdx >= numAngles){
		return;
	}


	int numberOfPixels[2];
	numberOfPixels[0]=numberOfPixelsX;
	numberOfPixels[1]=numberOfPixelsY;

	double pixSizes[2];
	pixSizes[0]=pixSizesX;
	pixSizes[1]=pixSizesY;

	int sinogramDataIdx = detPixIdx*numAngles + angIdx; //Column-major order!
	sinogramDataIdx =     angIdx*numDetPixels + detPixIdx;

	double halfPhantomWidth = numberOfPixels[0]*pixSizes[0]/2;
	double halfPhantomHeight = numberOfPixels[1]*pixSizes[1]/2;

	//trigonometric functions of the angles
	double theta = angles [angIdx];
	double sinTheta = sin(theta); //TODO: hasznaljuk a sincos fuggvenyt!!!!
	double cosTheta = cos(theta);

	//Distance of the detector plane from origin which is outside of the phantom
	double detDist = 1.1 * sqrt(pow(pixSizes[0]*numberOfPixels[0], 2) + pow(pixSizes[1]*numberOfPixels[1], 2) ); //Distance between the detector plane and centre of rotation

	const double* dataPtr=phantom;    //TODO: rename dataPtr -> phantom

	const double invPixSize0 = 1 / pixSizes[0];
	const double invPixSize1 = 1 / pixSizes[1];

	const double pixSizeRatio01 = pixSizes[0] / pixSizes[1];
	const double pixSizeRatio10 = pixSizes[1] / pixSizes[0];

	double p1[2];
	double p2[2];

	//beam intersects the columns at most two pixels

	if( pixSizes[1] / pixSizes[0] >= std::abs(std::tan(M_PI/2-theta)) ){

	    double t = -1*numDetPixels*detPixSize/2+(detPixIdx+0.5)*detPixSize;

	    p1[0]=detDist * sinTheta + t * cosTheta;
	    p1[1]=-1*detDist * cosTheta + t * sinTheta;

	    p2[0] = -1 * detDist * sinTheta + t * cosTheta;
	    p2[1] = detDist * cosTheta + t * sinTheta;

    	double ky = (p1[1]-p2[1])/(p1[0]-p2[0]);
    	double pathInSinglePixel = sqrt(1+ky*ky)*pixSizes[0];

    	//go through the columns of the image
    	for(int colIdx=0; colIdx<numberOfPixels[0]; ++colIdx){
    		double yi_minus = ( halfPhantomHeight - (ky*( colIdx   *pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) * invPixSize1; //Always on the left side of the column
    		//double yi_plus  = ( halfPhantomHeight - (ky*((colIdx+1)*pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) / pixSizes[1]; //Always on the right side of the column      //Optimized code in next line
    		double yi_plus = yi_minus - ky*pixSizeRatio01;

    		int Yi_minusIdx = floor(yi_minus);
	    	int Yi_plusIdx = floor(yi_plus);

	    	//int Yi_minusIdx = static_cast<int>(yi_minus) - ( yi_minus < static_cast<int>(yi_minus));  //Seemingly faster floor, but not
	    	//int Yi_plusIdx = static_cast<int>(yi_plus) - ( yi_plus < static_cast<int>(yi_plus));

	    	double l_minus, l_plus; //intersecting lengths when crossing two pixels
	    	if( Yi_minusIdx == Yi_plusIdx ){ //intersecting only one pixel
	    		if( (Yi_minusIdx < numberOfPixels[1]) and (Yi_minusIdx >= 0 ) ){
	    			//l=sqrt(1+ky*ky)*pixSizes[0]; //Optimized away with pathInSinglePixel

	    			sinogram[sinogramDataIdx] += pathInSinglePixel * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
	    		}
	    	}
	    	else{
	    		if ( (Yi_minusIdx < numberOfPixels[1]) and (Yi_minusIdx >= 0) ){
	    			l_minus=(max(Yi_minusIdx, Yi_plusIdx)-yi_minus) / (yi_plus - yi_minus) * pathInSinglePixel;

	    			sinogram[sinogramDataIdx] += l_minus * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
	    			}

	    		if ( (Yi_plusIdx < numberOfPixels[1]) and (Yi_plusIdx >= 0 ) ){
	    			l_plus=(yi_plus - max(Yi_minusIdx, Yi_plusIdx)) / (yi_plus - yi_minus) * pathInSinglePixel;

	    			sinogram[sinogramDataIdx] += l_plus * dataPtr[Yi_plusIdx*numberOfPixels[0] + colIdx];
	    		}
	    	}
	    }
	}
	else{      //beam intersects the rows at most two pixels

		double t = -1*numDetPixels*detPixSize/2+(detPixIdx+0.5)*detPixSize;

	    p1[0]=detDist * sinTheta + t * cosTheta;
	    p1[1]=-1*detDist * cosTheta + t * sinTheta;

	    p2[0] = -1 * detDist * sinTheta + t * cosTheta;
	    p2[1] = detDist * cosTheta + t * sinTheta;

	   	double kx = (p1[0]-p2[0])/(p1[1]-p2[1]);
	   	double pathInSinglePixel = sqrt(1+kx*kx)*pixSizes[1];

	   	//go through the rows of the image
	    for(int rowIdx=0; rowIdx<numberOfPixels[1]; ++rowIdx){

	    	double xi_minus = (halfPhantomWidth + (kx*( halfPhantomHeight - rowIdx     *pixSizes[1] - p1[1] ) + p1[0] ) )  * invPixSize0;
	    	//double xi_plus  = (halfPhantomWidth + (kx*( halfPhantomHeight - (rowIdx+1) *pixSizes[1] - p1[1] ) + p1[0] ) ) / pixSizes[0];    //Optimized code in next line
	    	double xi_plus = xi_minus - kx*pixSizeRatio10;

	    	int Xi_minusIdx = floor(xi_minus);
	        int Xi_plusIdx = floor(xi_plus);

	        //int Xi_minusIdx = static_cast<int>(xi_minus) - ( xi_minus < static_cast<int>(xi_minus));  //seemingly faster floor, but NOT
	        //int Xi_plusIdx = static_cast<int>(xi_plus) - ( xi_plus < static_cast<int>(xi_plus));

	        double l_minus, l_plus; //intersecting lengths
	        if( Xi_minusIdx == Xi_plusIdx ){
	        	if( (Xi_minusIdx < numberOfPixels[0]) and (Xi_minusIdx >= 0 ) ){
	        		//l=sqrt(1+kx*kx)*pixSizes[1]; //Optimized away with pathInSinglePixel

	        		sinogram[sinogramDataIdx] += pathInSinglePixel * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
	        	}
	        }
	        else{
	        	if ( (Xi_minusIdx < numberOfPixels[0]) and (Xi_minusIdx >= 0 ) ){
	        		l_minus=((Xi_minusIdx, Xi_plusIdx)-xi_minus) / (xi_plus - xi_minus) * pathInSinglePixel;

	        		sinogram[sinogramDataIdx] += l_minus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
	        	}

	    		if ( (Xi_plusIdx < numberOfPixels[0]) and (Xi_plusIdx >= 0) ){
	    			l_plus=(xi_plus - max(Xi_minusIdx, Xi_plusIdx)) / (xi_plus - xi_minus) * pathInSinglePixel;

	    			sinogram[sinogramDataIdx] += l_plus * dataPtr[rowIdx*numberOfPixels[0] + Xi_plusIdx];
	        	}
	        }
	    }
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
/*
	//Allocate space and copy data
	double *d_phantom;
	checkCudaErrors(cudaMalloc(&d_phantom, sizeof(double) * numberOfPixels[0]*numberOfPixels[1] ));
	checkCudaErrors(cudaMemcpy(d_phantom, actualPhantom.getDataAsEigenMatrixRef().data(),
			                   sizeof(double) * numberOfPixels[0]*numberOfPixels[1], cudaMemcpyHostToDevice));

	double *d_sinogram;
	checkCudaErrors(cudaMalloc(&d_sinogram, sizeof(double) * numAngles * pixNum ));
	checkCudaErrors(cudaMemset(d_sinogram, 0, sizeof(double) * numAngles * pixNum )); //TODO: Lehet h jobb/gyorsabb ha a kernelben nullazzuk?

	double *d_angles;
	checkCudaErrors(cudaMalloc(&d_angles, sizeof(double) * numAngles ));
	checkCudaErrors(cudaMemcpy(d_angles, angles.data(),
				                   sizeof(double) * numAngles, cudaMemcpyHostToDevice));

	//CALL THE PROJECTION KERNEL!!!!
	const dim3 blockSize(16,16);
	const dim3 gridSize(pixNum/blockSize.x+1,numAngles/blockSize.y+1);
	rayDrivenKernel<<<gridSize, blockSize>>>(d_phantom, numberOfPixels[0], numberOfPixels[1], pixSizes[0], pixSizes[1],
			                                 d_sinogram, numAngles, d_angles, pixNum, detWidth/pixNum);
	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	//Read back the result and free memory
	checkCudaErrors(cudaMemcpy(sinogram.data(), d_sinogram,
					                   sizeof(double) * numAngles*pixNum, cudaMemcpyDeviceToHost));

	cudaFree(d_sinogram);
	cudaFree(d_phantom);
	cudaFree(d_angles);
*/
	////////////////////////////////////////////////////////////////////////////////
	/// Projection with ray-driven method on GPU ENDS here !!
	////////////////////////////////////////////////////////////////////////////////

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
	std::cout << "Projection with ray-driven method took " << duration.count() << " milliseconds" << std::endl;

	return sinogram;
}

#endif
