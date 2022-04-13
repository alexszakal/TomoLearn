#include <Gen1CT.hpp>

#include <config.h>
#include <cstdio>
#include <chrono>

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
void rayDrivenKernel(const double* phantom, int numberOfPixelsX, int numberOfPixelsY, double pixSizesX, double pixSizesY,
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

	int sinogramDataIdx = angIdx*numDetPixels + detPixIdx; //Column-major order!

	double halfPhantomWidth = numberOfPixelsX*pixSizesX/2;
	double halfPhantomHeight = numberOfPixelsY*pixSizesY/2;

	//trigonometric functions of the angles
	double theta = angles [angIdx];
	double sinTheta, cosTheta;
	sincos(theta, &sinTheta, &cosTheta);

	//Distance of the detector plane from origin which is outside of the phantom
	double detDist = 1.1 * sqrt(pow(pixSizesX*numberOfPixelsX, 2) + pow(pixSizesY*numberOfPixelsY, 2) ); //Distance between the detector plane and centre of rotation

	double p1[2];
	double p2[2];

	double sinoPointValue=0.0;  //Temporary local variable to accumulate result

	//beam intersects the columns at most two pixels
	if( pixSizesY / pixSizesX >= std::abs(std::tan(M_PI/2-theta)) ){
		const double invPixSize1 = 1 / pixSizesY;
		const double pixSizeRatio01 = pixSizesX / pixSizesY;

	    const double t = -1*numDetPixels*detPixSize/2+(detPixIdx+0.5)*detPixSize;

	    p1[0]=detDist * sinTheta + t * cosTheta;
	    p1[1]=-1*detDist * cosTheta + t * sinTheta;

	    p2[0] = -1 * detDist * sinTheta + t * cosTheta;
	    p2[1] = detDist * cosTheta + t * sinTheta;

    	const double ky = (p1[1]-p2[1])/(p1[0]-p2[0]);
    	const double pathInSinglePixel = sqrt(1+ky*ky)*pixSizesX;

    	double yi_minus = ( halfPhantomHeight - (ky*( (-1)   *pixSizesX - halfPhantomWidth - p1[0] ) + p1[1] ) ) * invPixSize1;
    	const double yi_minusIncrement = -1*ky*invPixSize1*pixSizesX;
    	const double yi_plusIncrement = -1 * ky*pixSizeRatio01;
    	//go through the columns of the image
    	for(int colIdx=0; colIdx<numberOfPixelsX; ++colIdx){
    		//double yi_minus = ( halfPhantomHeight - (ky*( colIdx   *pixSizesX - halfPhantomWidth - p1[0] ) + p1[1] ) ) * invPixSize1; //Always on the left side of the column
    		//double yi_plus  = ( halfPhantomHeight - (ky*((colIdx+1)*pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) / pixSizes[1]; //Always on the right side of the column      //Optimized code in next line
    		//double yi_plus = yi_minus - ky*pixSizeRatio01;
    		yi_minus += yi_minusIncrement;
    		double yi_plus = yi_minus + yi_plusIncrement;

    		int Yi_minusIdx = floor(yi_minus);
	    	int Yi_plusIdx = floor(yi_plus);

	    	//int Yi_minusIdx = static_cast<int>(yi_minus) - ( yi_minus < static_cast<int>(yi_minus));  //Seemingly faster floor, but not
	    	//int Yi_plusIdx = static_cast<int>(yi_plus) - ( yi_plus < static_cast<int>(yi_plus));

	    	double l_minus, l_plus; //intersecting lengths when crossing two pixels
	    	if( Yi_minusIdx == Yi_plusIdx ){ //intersecting only one pixel
	    		if( (Yi_minusIdx < numberOfPixelsY) and (Yi_minusIdx >= 0 ) ){
	    			//l=sqrt(1+ky*ky)*pixSizes[0]; //Optimized away with pathInSinglePixel

	    			sinoPointValue += pathInSinglePixel * phantom[Yi_minusIdx*numberOfPixelsX + colIdx];
	    		}
	    	}
	    	else{
	    		if ( (Yi_minusIdx < numberOfPixelsY) and (Yi_minusIdx >= 0) ){
	    			l_minus=(max(Yi_minusIdx, Yi_plusIdx)-yi_minus) / (yi_plus - yi_minus) * pathInSinglePixel;

	    			sinoPointValue += l_minus * phantom[Yi_minusIdx*numberOfPixelsX + colIdx];
	    			//We have l_minus -> we can calculate l_plus with only a subtraction from pathInSinglePixel
	    			if( (Yi_plusIdx < numberOfPixelsY) and (Yi_plusIdx >= 0 ) ){
	    				sinoPointValue += (pathInSinglePixel- l_minus) * phantom[Yi_plusIdx*numberOfPixelsX + colIdx];
	    				continue;
	    			}
	    		}

	    		if ( (Yi_plusIdx < numberOfPixelsY) and (Yi_plusIdx >= 0 ) ){
	    			l_plus=(yi_plus - max(Yi_minusIdx, Yi_plusIdx)) / (yi_plus - yi_minus) * pathInSinglePixel;

	    			sinoPointValue += l_plus * phantom[Yi_plusIdx*numberOfPixelsX + colIdx];
	    		}
	    	}
	    }
	}
	else{      //beam intersects the rows at most two pixels
		const double invPixSize0 = 1 / pixSizesX;
		const double pixSizeRatio10 = pixSizesY / pixSizesX;

		const double t = -1*numDetPixels*detPixSize/2+(detPixIdx+0.5)*detPixSize;

	    p1[0]=detDist * sinTheta + t * cosTheta;
	    p1[1]=-1*detDist * cosTheta + t * sinTheta;

	    p2[0] = -1 * detDist * sinTheta + t * cosTheta;
	    p2[1] = detDist * cosTheta + t * sinTheta;

	   	const double kx = (p1[0]-p2[0])/(p1[1]-p2[1]);
	   	const double pathInSinglePixel = sqrt(1+kx*kx)*pixSizesY;

	   	//go through the rows of the image
	   	double xi_minus = (halfPhantomWidth + (kx*( halfPhantomHeight - (-1)*pixSizesY - p1[1] ) + p1[0] ) )  * invPixSize0;
	   	const double xi_minusIncrement = -1*kx*invPixSize0*pixSizesY;
	   	const double xi_plusIncrement = -1*kx*pixSizeRatio10;
	    for(int rowIdx=0; rowIdx<numberOfPixelsY; ++rowIdx){
	    	//double xi_plus  = (halfPhantomWidth + (kx*( halfPhantomHeight - (rowIdx+1) *pixSizes[1] - p1[1] ) + p1[0] ) ) / pixSizes[0];    //Optimized code in next line
	    	//double xi_plus = xi_minus - kx*pixSizeRatio10;
	    	xi_minus += xi_minusIncrement;
	    	double xi_plus = xi_minus + xi_plusIncrement;

	    	int Xi_minusIdx = floor(xi_minus);
	        int Xi_plusIdx = floor(xi_plus);

	        //int Xi_minusIdx = static_cast<int>(xi_minus) - ( xi_minus < static_cast<int>(xi_minus));  //seemingly faster floor, but NOT
	        //int Xi_plusIdx = static_cast<int>(xi_plus) - ( xi_plus < static_cast<int>(xi_plus));

	        double l_minus, l_plus; //intersecting lengths
	        if( Xi_minusIdx == Xi_plusIdx ){
	        	if( (Xi_minusIdx < numberOfPixelsX) and (Xi_minusIdx >= 0 ) ){
	        		//l=sqrt(1+kx*kx)*pixSizes[1]; //Optimized away with pathInSinglePixel

	        		sinoPointValue += pathInSinglePixel * phantom[rowIdx*numberOfPixelsX + Xi_minusIdx];
	        	}
	        }
	        else{
	        	if ( (Xi_minusIdx < numberOfPixelsX) and (Xi_minusIdx >= 0 ) ){
	        		l_minus=(max(Xi_minusIdx, Xi_plusIdx)-xi_minus) / (xi_plus - xi_minus) * pathInSinglePixel;

	        		sinoPointValue += l_minus * phantom[rowIdx*numberOfPixelsX + Xi_minusIdx];
	        		//We have l_minus -> we can calculate l_plus with only a subtraction from pathInSinglePixel
	        		if ( (Xi_plusIdx <= numberOfPixelsX) and (Xi_plusIdx >= 0 ) ){   //Shortcut to avoid one more std::max call
	        			sinoPointValue += (pathInSinglePixel-l_minus) * phantom[rowIdx*numberOfPixelsX + Xi_minusIdx];
	        			continue;
	        		}
	        	}

	    		if ( (Xi_plusIdx <= numberOfPixelsX) and (Xi_plusIdx >= 0) ){
	    			l_plus=(xi_plus - max(Xi_minusIdx, Xi_plusIdx)) / (xi_plus - xi_minus) * pathInSinglePixel;

	    			sinoPointValue += l_plus * phantom[rowIdx*numberOfPixelsX + Xi_plusIdx];
	        	}
	        }
	    }
	}

	sinogram[sinogramDataIdx] = sinoPointValue; //Write the result to global memory

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

void launchRayDrivenKernel(const double* phantomData, std::array<int, 2> numberOfPixels, std::array<double, 2> pixSizes,
		                   int numAngles, const double* anglesData, int pixNum, double detWidth,
						   double* sinogramData){

	auto start = std::chrono::high_resolution_clock::now();

	//Allocate space and copy data
	double *d_phantom;
	checkCudaErrors(cudaMalloc(&d_phantom, sizeof(double) * numberOfPixels[0]*numberOfPixels[1] ));
	checkCudaErrors(cudaMemcpy(d_phantom, phantomData,
			                   sizeof(double) * numberOfPixels[0]*numberOfPixels[1], cudaMemcpyHostToDevice));

	double *d_sinogram;
	checkCudaErrors(cudaMalloc(&d_sinogram, sizeof(double) * numAngles * pixNum ));

	double *d_angles;
	checkCudaErrors(cudaMalloc(&d_angles, sizeof(double) * numAngles ));
	checkCudaErrors(cudaMemcpy(d_angles, anglesData,
					                   sizeof(double) * numAngles, cudaMemcpyHostToDevice));

	auto phase1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(phase1 - start);
	std::cout << "\nMemory allocation and data transfer to GPU took " << duration.count() << " milliseconds";

	//CALL THE PROJECTION KERNEL!!!!
	const dim3 blockSize(16,16);
	const dim3 gridSize(pixNum/blockSize.x+1,numAngles/blockSize.y+1);
	rayDrivenKernel<<<gridSize, blockSize>>>(d_phantom, numberOfPixels[0], numberOfPixels[1], pixSizes[0], pixSizes[1],
				                                 d_sinogram, numAngles, d_angles, pixNum, detWidth/pixNum);
	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	auto phase2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds >(phase2 - phase1);
	std::cout << "\n Running the rayDriven projection kernel on GPU took " << duration.count() << " milliseconds";

	//Read back the result and free memory
	checkCudaErrors(cudaMemcpy(sinogramData, d_sinogram,
					                   sizeof(double) * numAngles*pixNum, cudaMemcpyDeviceToHost));

	cudaFree(d_sinogram);
	cudaFree(d_phantom);
	cudaFree(d_angles);

	auto phase3 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds >(phase3 - phase2);
	std::cout << "\nCopy back results and free GPU memory took " << duration.count() << " milliseconds";

}

#endif
