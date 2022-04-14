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
void rayDrivenProjectionKernel(const double* phantom, int numberOfPixelsX, int numberOfPixelsY, double pixSizesX, double pixSizesY,
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
    		yi_minus += yi_minusIncrement;
    		double yi_plus = yi_minus + yi_plusIncrement;

    		int Yi_minusIdx = floor(yi_minus);
	    	int Yi_plusIdx = floor(yi_plus);

	    	double l_minus, l_plus; //intersecting lengths when crossing two pixels
	    	if( Yi_minusIdx == Yi_plusIdx ){ //intersecting only one pixel
	    		if( (Yi_minusIdx < numberOfPixelsY) and (Yi_minusIdx >= 0 ) ){
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
	    	xi_minus += xi_minusIncrement;
	    	double xi_plus = xi_minus + xi_plusIncrement;

	    	int Xi_minusIdx = floor(xi_minus);
	        int Xi_plusIdx = floor(xi_plus);

	        double l_minus, l_plus; //intersecting lengths
	        if( Xi_minusIdx == Xi_plusIdx ){
	        	if( (Xi_minusIdx < numberOfPixelsX) and (Xi_minusIdx >= 0 ) ){
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

	    		if ( (Xi_plusIdx < numberOfPixelsX) and (Xi_plusIdx >= 0) ){
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

void launchRayDrivenProjectionKernel(const double* phantomData, std::array<int, 2> numberOfPixels, std::array<double, 2> pixSizes,
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
	rayDrivenProjectionKernel<<<gridSize, blockSize>>>(d_phantom, numberOfPixels[0], numberOfPixels[1], pixSizes[0], pixSizes[1],
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

__global__
void rayDrivenBackprojectionKernel(double* d_sinogram, int numAngles, double* d_angles, double* sinThetaVector, double* cosThetaVector, int pixNum, double detWidth,
									double* d_backProjection, int numberOfPixelsX, int numberOfPixelsY, double resolutionX, double resolutionY){

	int xIdx = blockIdx.x*blockDim.x + threadIdx.x;
	int yIdx = blockIdx.y*blockDim.y + threadIdx.y;

	if ( xIdx >= numberOfPixelsX){
		return;
	}
	if ( yIdx >= numberOfPixelsY){
		return;
	}

	int backProjectionDataIdx = yIdx*numberOfPixelsX + xIdx; //Column-major order!

	const double pixelRadius = sqrt(pow(resolutionX/2,2) + pow(resolutionY/2,2)); //Radius of circle drawn around a pixel
    const double halfDetWidth = detWidth / 2;
    const double invDetPixSize= pixNum/detWidth;

    double xValue= -1*numberOfPixelsX*resolutionX/2 + resolutionX/2 + xIdx*resolutionX;
    double yValue=    numberOfPixelsY*resolutionY/2 - resolutionY/2 - yIdx*resolutionY;

    const double pixelSize = detWidth/pixNum;

    d_backProjection[backProjectionDataIdx] = 0.0;

	//Calculate the value in pixel
	//Go through the angles
	for(int angIdx=0; angIdx<numAngles; ++angIdx){
		//Determine the contributing detector pixels (i.e. rays)
		double xr = xValue*cosThetaVector[static_cast<size_t>(angIdx)]+
				       yValue*sinThetaVector[static_cast<size_t>(angIdx)];   //This corresponds to "t" of the image pixel center

		double minPixIdx = (xr + halfDetWidth - pixelRadius )*invDetPixSize;
		double maxPixIdx = (xr + halfDetWidth + pixelRadius )*invDetPixSize;

		//Calculate the intersection length and add to the backprojected image
		//Go through the possible pixels
		double lSum=0;
		double angleContrib = 0;
		double a = cosThetaVector[static_cast<size_t>(angIdx)];
		double b = sinThetaVector[static_cast<size_t>(angIdx)];
		double absa = abs(a);
		double absb = abs(b);

		//Calculate the boundaries of pixel intersection types
		double d_max = (std::abs(a*resolutionX) + std::abs(b*resolutionY))/2; //For definition see HaoGao's article
		double d_min =  std::abs( std::abs(a*resolutionX) - std::abs(b*resolutionY) )/2;

		for(int detPixIdx=max(0,static_cast<int>(minPixIdx));
			    detPixIdx <= min(static_cast<int>(maxPixIdx)+0, static_cast<int>(pixNum-1) );
			    ++detPixIdx){
			double c = -1*halfDetWidth+(detPixIdx+0.5)*pixelSize;  // "t" of the detector pixel center

			//Calculate the intersection length
			double d_act=std::abs(a*xValue + b*yValue - c);
			double l;

			if(d_act < d_min){
				if( absa < absb ){
					l=resolutionX/absb;
				}
				else{
					l=resolutionY/absa;
				}
			}
			else if( (d_act >= d_min) and (d_act < d_max) ){
				l=(d_max-d_act)/absa/absb;
			}
			else{
				l=0;
			}

			lSum += l;
			angleContrib += d_sinogram[angIdx*pixNum+detPixIdx]*l;
		}
		if(lSum != 0){
			d_backProjection[backProjectionDataIdx] += angleContrib/lSum;
		}
	}

	//Multiply with dTheta
	d_backProjection[backProjectionDataIdx] = d_backProjection[backProjectionDataIdx]*M_PI/numAngles;
}

void launchRayDrivenBackprojectionKernel(const double* sinogram, int numAngles, const double* anglesData, int pixNum, double detWidth,
										 double* backProjection, const std::array<int,2>& numberOfPixels, const std::array<double,2>& resolution){

	auto start = std::chrono::high_resolution_clock::now();

	double *d_sinogram;
	checkCudaErrors(cudaMalloc(&d_sinogram, sizeof(double) * numAngles*pixNum ));
	checkCudaErrors(cudaMemcpy(d_sinogram, sinogram,
				                   sizeof(double) * numAngles*pixNum, cudaMemcpyHostToDevice));

	double *d_angles;
	checkCudaErrors(cudaMalloc(&d_angles, sizeof(double) * numAngles ));
	checkCudaErrors(cudaMemcpy(d_angles, anglesData,
					                   sizeof(double) * numAngles, cudaMemcpyHostToDevice));

	double *d_backprojection;
	checkCudaErrors(cudaMalloc(&d_backprojection, sizeof(double) * numberOfPixels[0]*numberOfPixels[1] ));

	std::vector<double>	thetaVector,
	                    sinThetaVector,
						cosThetaVector;

	//Vector of theta values and trigonometric functions
	for(int i=0; i<numAngles; i++){
	  	thetaVector.push_back( std::fmod(anglesData[i], 2*M_PI) );
	   	sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
	   	cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
	}

	double *d_sinThetaVector;
	checkCudaErrors(cudaMalloc(&d_sinThetaVector, sizeof(double) * numAngles ));
	checkCudaErrors(cudaMemcpy(d_sinThetaVector, sinThetaVector.data(),
						                   sizeof(double) * numAngles, cudaMemcpyHostToDevice));

	double *d_cosThetaVector;
	checkCudaErrors(cudaMalloc(&d_cosThetaVector, sizeof(double) * numAngles ));
	checkCudaErrors(cudaMemcpy(d_cosThetaVector, cosThetaVector.data(),
						                   sizeof(double) * numAngles, cudaMemcpyHostToDevice));

	auto stop1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start);
	std::cout << "Phase1 took " << duration.count() << " milliseconds" << std::endl;

	auto start2 = std::chrono::high_resolution_clock::now();
	const dim3 blockSize(16,16);
	const dim3 gridSize(numberOfPixels[0]/blockSize.x+1,numberOfPixels[1]/blockSize.y+1);
	rayDrivenBackprojectionKernel<<<gridSize, blockSize>>>(d_sinogram, numAngles, d_angles, d_sinThetaVector, d_cosThetaVector, pixNum, detWidth,
			d_backprojection, numberOfPixels[0], numberOfPixels[1], resolution[0], resolution[1]);
	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	auto stop2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
	std::cout << "Running the kernel took " << duration.count() << " milliseconds" << std::endl;

	auto start3 = std::chrono::high_resolution_clock::now();
	//Read back the result and free memory
	checkCudaErrors(cudaMemcpy(backProjection, d_backprojection,
						         sizeof(double) * numberOfPixels[0]*numberOfPixels[1], cudaMemcpyDeviceToHost));

	cudaFree(d_sinogram);
	cudaFree(d_angles);
	cudaFree(d_backprojection);

	auto stop3 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop3 - start3);
	std::cout << "Reading the results and free memory took " << duration.count() << " milliseconds" << std::endl;

}
#endif




