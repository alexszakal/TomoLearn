
#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

#include <Gen1CT.hpp>
#include <Phantom.hpp>
#include <Reconst.hpp>
#include <CTScan.hpp>
#include <Filter.hpp>
#include <cudaFunctions.cuh>

#include <matplot/matplot.h>

#include <utility>   //std::swap
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <chrono>
#include <array>
#include <vector>
#include <random>
#include <algorithm> //std::min and std::max and std::clamp

#include <config.h>


/***
 * Default constructor; width of detector: 100mm, number of pixels: 100
 */
Gen1CT::Gen1CT():detWidth{100},pixNum{100}{
	pixPositions.resize( static_cast<size_t>(pixNum) );
	double pixelSize{detWidth/pixNum};
	for (int i=0; i<pixNum; i++){
		pixPositions[ static_cast<size_t>(i) ]=-1*detWidth/2+(i+0.5)*pixelSize;
	}
};

/***
 *Gen1CT constructor
 * @param detWidth Width of the detector [mm]
 * @param pixNum Number of detector pixels
 */
Gen1CT::Gen1CT(double detWidth, int pixNum):detWidth{detWidth},pixNum{pixNum}{
	pixPositions.resize( static_cast<size_t>(pixNum) );
	double pixelSize{detWidth/pixNum};
	for (int i=0; i<pixNum; i++){
		pixPositions[ static_cast<size_t>(i) ]=-1*detWidth/2+(i+0.5)*pixelSize;
	}
};

/**
 * Add a phantom from an existing image
 * @param label The label of the phantom in the library
 * @param phantomImageSource Location of the image file
 * @param pixSizes Pixel size in [mm]
 * @param convertFromHUtoLA The values in phantomImagesource are in Hounsfield units, conversion is needed to linear attenuation units.
 */
void Gen1CT::addPhantom(const std::string &label,
		const std::string &phantomImageSource,
		std::array<double, 2> pixSizes,
		bool convertFromHUtoLA) {

	auto it = phantoms.find(label);
	if (it != phantoms.end()) {
		std::cout << std::endl << "WARNING! A phantom with label \"" << label
				<< "\"is already loaded!!! Overwriting!!!";
		phantoms.erase(it);
	}
	phantoms.emplace(label, Phantom(label, phantomImageSource, pixSizes, convertFromHUtoLA));
}

/**
 * Add an existing Phantom object to the Phantom library of the Gen1CT object
 * @param newPhantom Phantom object to add to the library
 */
void Gen1CT::addPhantom(const Phantom &newPhantom) {
	std::string phantomLabel = newPhantom.getLabel();
	auto it = phantoms.find(phantomLabel);
	if (it != phantoms.end()) {
		std::cout << std::endl << "WARNING! A phantom with label \""
				<< phantomLabel << "\"is already loaded!!! Overwriting!!!";
		phantoms.erase(it);
	}
	phantoms.emplace(phantomLabel, Phantom(newPhantom));
}

/***
 * Interactive display of the Phantom with 'label' from the Phantoms library
 * @param label Label of the Phantom to be displayed
 */
void Gen1CT::displayPhantom(const std::string& label){
	if(phantoms.find(label) != phantoms.end()){
		phantoms[label].display(label);
	}
	else
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

/***
 * Setter function for the I0 variable which sets the number of photons starting from the tube at a measurement point
 * @param newI0 The new value of I0
 */
void Gen1CT::setI0(double newI0){
	I0=newI0;
}

/***
 * Simulate the measurement of a Phantom
 * @param phantomLabel Label of the Phantom in the Phantom library
 * @param angles Angles of the measurement [rad]
 * @param scanLabel Label of the result in the library for the scans
 * @param projector Type of the applied projector
 */
void Gen1CT::measure(const std::string& phantomLabel,
		     const Eigen::VectorXd& angles,
			 const std::string& scanLabel,
			 projectorType projector){

	if(phantoms.find(phantomLabel) == phantoms.end()){
		std::cout << std::endl << "ERROR!! phantomLabel: \"" << phantomLabel << "\" could not be found!! Abort mission";
		return;
	}

	Phantom& actualPhantom = phantoms[phantomLabel];

	Eigen::MatrixXd sinogram = project(actualPhantom, angles, projector);

	//Check if scanLabel can be found in CTScans map
	auto it = scans.find(scanLabel);
	if(it != scans.end()){
		std::cout << std::endl << "WARNING! A scan with label \"" << scanLabel << "\" already exists!!! Overwriting!!!";
		scans.erase(it);
	}

	if(I0 !=0.0){ //Simulation of the effect of statistics on the line integrals has to be done
		//Calculate the expected value
		sinogram = Eigen::exp( sinogram.array()* (-1.0) ) * I0 ;

		//Randomize the matrix
		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		std::mt19937_64 generator(seed1);    //Seed the generator

		for(int i=0; i<sinogram.rows(); i++){
			for(int j=0; j<sinogram.cols(); j++){
				std::poisson_distribution<int> poissonDist( sinogram(i,j) );
				sinogram(i,j) = poissonDist(generator);
			}
		}
	}

	//Save the scan in the library
	scans.emplace(scanLabel, CTScan(scanLabel,sinogram, detWidth, angles, I0));
}

/***
 * Calculate the projection of a Phantom object
 * @param actualPhantom The Phantom object of which the projetion is calculated
 * @param angles Projection angles [rad]
 * @param projector Type of the used projector
 * @return The sinogram matrix of the actualPhantom
 */
Eigen::MatrixXd Gen1CT::project(const Phantom& actualPhantom, const Eigen::VectorXd& angles, projectorType projector){
	Eigen::MatrixXd sinogram;

	switch (projector){
		case projectorType::pixelDriven:
			sinogram = project_pixelDriven_CPU(actualPhantom, angles);
			break;
		case projectorType::Siddon:
			sinogram = project_Siddon_CPU(actualPhantom, angles);
			break;
		case projectorType::rayDriven:
			sinogram = project_rayDriven_CPU(actualPhantom, angles);
			break;
#if ENABLE_CUDA
		case projectorType::rayDriven_GPU:
			sinogram = project_rayDriven_GPU(actualPhantom, angles);
			break;
#endif
		case projectorType::rayDrivenOptimized:
			sinogram = project_rayDrivenOptimized_CPU(actualPhantom, angles);
			break;
	}

	return sinogram;
}

/***
 * Ray-driven projection calculated on the CPU
 * @param actualPhantom The Phantom object which has to be projected
 * @param angles Projection angles [rad]
 * @return The sinogram matrix of the actualPhantom
 */
Eigen::MatrixXd Gen1CT::project_rayDriven_CPU(const Phantom& actualPhantom,
		                    const Eigen::VectorXd& angles){
	std::cout <<  "\nProjection with ray-driven method started";
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(static_cast<long>(pixNum), static_cast<long>(numAngles));

	auto pixSizes = actualPhantom.getPixSizes();
    auto numberOfPixels = actualPhantom.getNumberOfPixels();

    ////////////////////////////////////////////////////////////////////////////////
    /// Projection with ray-driven method STARTS here !!
    ////////////////////////////////////////////////////////////////////////////////

    double halfPhantomWidth = numberOfPixels[0]*pixSizes[0]/2;
    double halfPhantomHeight = numberOfPixels[1]*pixSizes[1]/2;

    //trigonometric functions of the angles
    std::vector<double>	thetaVector,
    	                sinThetaVector,
						cosThetaVector;

    for(int i=0; i<numAngles; i++){
    	thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
    	sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
    	cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
    }

    //Distance of the detector plane from origin which is outside of the phantom
    double detDist = 1.1 * std::sqrt(std::pow(pixSizes[0]*numberOfPixels[0], 2) + pow(pixSizes[1]*numberOfPixels[1], 2) ); //Distance between the detector plane and centre of rotation

    const double* dataPtr=actualPhantom.getDataAsEigenMatrixRef().data();

    const double invPixSize0 = 1 / pixSizes[0];
    const double invPixSize1 = 1 / pixSizes[1];

    const double pixSizeRatio01 = pixSizes[0] / pixSizes[1];
    const double pixSizeRatio10 = pixSizes[1] / pixSizes[0];

    double p1[2];
    double p2[2];
    //Go through the angles
    for(int angI=0; angI < numAngles; ++angI){
    	//beam intersects the columns at most two pixels
    	if( pixSizes[1] / pixSizes[0] >= std::abs(std::tan(M_PI/2-thetaVector[static_cast<size_t>(angI)])) ){
    		//go through the different t values
    		for(int pixIdx=0; pixIdx<pixNum; ++pixIdx){
    			double t = pixPositions[static_cast<size_t>(pixIdx)];

    			p1[0]=detDist * sinThetaVector[static_cast<size_t>(angI)] + t * cosThetaVector[static_cast<size_t>(angI)];
    			p1[1]=-1*detDist * cosThetaVector[static_cast<size_t>(angI)] + t * sinThetaVector[static_cast<size_t>(angI)];

    			p2[0] = -1 * detDist * sinThetaVector[static_cast<size_t>(angI)] + t * cosThetaVector[static_cast<size_t>(angI)];
    			p2[1] = detDist * cosThetaVector[static_cast<size_t>(angI)] + t * sinThetaVector[static_cast<size_t>(angI)];

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
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += pathInSinglePixel * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
    					}
    				} else{
    					if ( (Yi_minusIdx < numberOfPixels[1]) and (Yi_minusIdx >= 0) ){
    						l_minus=(std::max(Yi_minusIdx, Yi_plusIdx)-yi_minus) / (yi_plus - yi_minus) * pathInSinglePixel;
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_minus * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
    						}

    					if ( (Yi_plusIdx < numberOfPixels[1]) and (Yi_plusIdx >= 0 ) ){
    						l_plus=(yi_plus - std::max(Yi_minusIdx, Yi_plusIdx)) / (yi_plus - yi_minus) * pathInSinglePixel;
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_plus * dataPtr[Yi_plusIdx*numberOfPixels[0] + colIdx];
    					}
    				}
    			}
    		}
    	}
    	else{      //beam intersects the rows at most two pixels
    	    //go through the different t values
    	    for(int pixIdx=0; pixIdx<pixNum; ++pixIdx){
    	    	double t = pixPositions[static_cast<size_t>(pixIdx)];

    	    	p1[0]=detDist * sinThetaVector[static_cast<size_t>(angI)] + t * cosThetaVector[static_cast<size_t>(angI)];
    	    	p1[1]=-1*detDist * cosThetaVector[static_cast<size_t>(angI)] + t * sinThetaVector[static_cast<size_t>(angI)];

    	    	p2[0] = -1 * detDist * sinThetaVector[static_cast<size_t>(angI)] + t * cosThetaVector[static_cast<size_t>(angI)];
    	    	p2[1] = detDist * cosThetaVector[static_cast<size_t>(angI)] + t * sinThetaVector[static_cast<size_t>(angI)];

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
    	    				sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += pathInSinglePixel * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
    	    			}
    	    		}
    	    		else{
    	    			if ( (Xi_minusIdx < numberOfPixels[0]) and (Xi_minusIdx >= 0 ) ){
    	    				l_minus=(std::max(Xi_minusIdx, Xi_plusIdx)-xi_minus) / (xi_plus - xi_minus) * pathInSinglePixel;
    	    				sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_minus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
    	    			}

    					if ( (Xi_plusIdx < numberOfPixels[0]) and (Xi_plusIdx >= 0) ){
    						l_plus=(xi_plus - std::max(Xi_minusIdx, Xi_plusIdx)) / (xi_plus - xi_minus) * pathInSinglePixel;
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_plus * dataPtr[rowIdx*numberOfPixels[0] + Xi_plusIdx];
    	    			}
    	    		}
    			}
    	    }
    	}
    }
    ////////////////////////////////////////////////////////////////////////////////
    /// Projection with ray-driven method ENDS here !!
    ////////////////////////////////////////////////////////////////////////////////

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedTime = stop - start;
    std::cout << "\nProjection with ray-driven method took " << elapsedTime.count() << " milliseconds\n";

    return sinogram;
}

/***
 * Ray-driven projection optimized for compute time calculated on the CPU
 * @param actualPhantom The Phantom object which has to be projected
 * @param angles Projection angles [rad]
 * @return The sinogram matrix of the actualPhantom
 */
Eigen::MatrixXd Gen1CT::project_rayDrivenOptimized_CPU(const Phantom& actualPhantom,
		                    const Eigen::VectorXd& angles){
	std::cout << "\nProjection with ray-driven method (optimized) started";
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(static_cast<long>(pixNum), static_cast<long>(numAngles));
	double* sinogramPtr = sinogram.data();

	auto pixSizes = actualPhantom.getPixSizes();
    auto numberOfPixels = actualPhantom.getNumberOfPixels();
    double invPixSize= pixNum / detWidth;

    ////////////////////////////////////////////////////////////////////////////////
    /// Projection with ray-driven method STARTS here !!
    ////////////////////////////////////////////////////////////////////////////////

    double halfPhantomWidth = numberOfPixels[0]*pixSizes[0]/2;
    double halfPhantomHeight = numberOfPixels[1]*pixSizes[1]/2;

    //trigonometric functions of the angles
    std::vector<double>	thetaVector,
    	                sinThetaVector,
						cosThetaVector,
						minusTanThetaVector,
    					minusCotanThetaVector;
    std::vector<bool>	intersectColumns;
    std::vector<double> pathInSinglePixelX,
	                    pathInSinglePixelY;

    for(int i=0; i<numAngles; i++){
    	thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
    	sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
    	cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
    	minusTanThetaVector.push_back( -1*tan(thetaVector[static_cast<size_t>(i)]));
    	minusCotanThetaVector.push_back(-1*1/tan(thetaVector[static_cast<size_t>(i)]));

    	pathInSinglePixelX.push_back(sqrt(1+minusCotanThetaVector[static_cast<size_t>(i)]*minusCotanThetaVector[static_cast<size_t>(i)])*pixSizes[0]);
    	pathInSinglePixelY.push_back(sqrt(1+minusTanThetaVector[static_cast<size_t>(i)]*minusTanThetaVector[static_cast<size_t>(i)])*pixSizes[1]);


    	intersectColumns.push_back( pixSizes[1] / pixSizes[0] >= std::abs(std::tan(M_PI/2-thetaVector[static_cast<size_t>(i)])) );
    }

    //Distance of the detector plane from origin which is outside of the phantom
    double detDist = 1.1 * sqrt(pow(pixSizes[0]*numberOfPixels[0], 2) + pow(pixSizes[1]*numberOfPixels[1], 2) );

    const double* dataPtr=actualPhantom.getDataAsEigenMatrixRef().data();

    const double invPixSize0 = 1 / pixSizes[0];
    const double invPixSize1 = 1 / pixSizes[1];

    const double pixSizeRatio01 = pixSizes[0] / pixSizes[1];
    const double pixSizeRatio10 = pixSizes[1] / pixSizes[0];

    double p1[2];
    //double p2[2]; //Not needed for the optimized calculation

    int subImageSize=128;  //This has to be tuned for the processor. 128 or 256 is ideal for the old ThinkPad
    double subImageRadius = std::sqrt( std::pow(subImageSize*pixSizes[0],2) + std::pow(subImageSize*pixSizes[1],2) );
    int subImageDimX=numberOfPixels[0]/subImageSize +1;
    int subImageDimY=numberOfPixels[1]/subImageSize +1;

    double leftmostPixCenter=(numberOfPixels[0]/2-0.5)*pixSizes[0];
    double bottomPixCenter=(numberOfPixels[1]/2-0.5)*pixSizes[1];
    //Go through the SubImages
    for(size_t subImageXIdx=0; subImageXIdx<static_cast<size_t>(subImageDimX); ++subImageXIdx){
    	for(size_t subImageYIdx=0; subImageYIdx<static_cast<size_t>(subImageDimY); ++subImageYIdx){

    		int colIdxMin=static_cast<int>(subImageXIdx)*subImageSize;
    		if(colIdxMin >= numberOfPixels[0])
    		    colIdxMin = numberOfPixels[0]-1;
    		int colIdxMax=colIdxMin+subImageSize-1;
    		if(colIdxMax >= numberOfPixels[0])
    			colIdxMax = numberOfPixels[0]-1;
    		int rowIdxMin=static_cast<int>(subImageYIdx)*subImageSize;
    		if(rowIdxMin >= numberOfPixels[1])
    		    rowIdxMin = numberOfPixels[1]-1;
    		int rowIdxMax=rowIdxMin+subImageSize-1;
    		if(rowIdxMax >= numberOfPixels[1])
    		    rowIdxMax = numberOfPixels[1]-1;

    		//Center of the subImage in laboratory coordinate system
    		double colMinPos = colIdxMin*pixSizes[0] - leftmostPixCenter;
    		double colMaxPos = colIdxMax*pixSizes[0] - leftmostPixCenter;
    		double subImageCenterX = (colMinPos + colMaxPos)/2 ;
    		double rowMinPos = bottomPixCenter - rowIdxMax*pixSizes[1];
    		double rowMaxPos = bottomPixCenter - rowIdxMin*pixSizes[1];
    		double subImageCenterY = (rowMinPos + rowMaxPos)/2 ;

    		//Go through the angles
			for(size_t angI=0; angI < static_cast<size_t>(numAngles); ++angI){
				//t-values that can intercept the subImage
				double subImageCenter_t = subImageCenterX * cosThetaVector[angI] + subImageCenterY * sinThetaVector[angI];
				int minPixIdx = std::ceil( (subImageCenter_t - subImageRadius + detWidth/2) * invPixSize );
				int maxPixIdx = minPixIdx + subImageRadius*2 * invPixSize ;

				//minPixIdx = std::clamp(minPixIdx, 0, static_cast<int>(pixNum-1));    //This is slower :(
				//maxPixIdx = std::clamp(maxPixIdx, 0, static_cast<int>(pixNum-1));
				if(minPixIdx<0)
					minPixIdx=0;
				if(minPixIdx>=pixNum)
					minPixIdx=pixNum-1;
				if(maxPixIdx<0)
					maxPixIdx=0;
				if(maxPixIdx>=pixNum)
					maxPixIdx=pixNum-1;

				//beam intersects the columns at most two pixels
				if( intersectColumns[angI] ){
					//go through the different t values
					for(size_t pixIdx=static_cast<size_t>(minPixIdx); pixIdx<=static_cast<size_t>(maxPixIdx); ++pixIdx){
						const double t = pixPositions[pixIdx];

						size_t sinogramPoint = pixIdx + angI*static_cast<size_t>(pixNum);       //Index of the point in sinogram array

						p1[0]=detDist * sinThetaVector[angI] + t * cosThetaVector[angI];
						p1[1]=-1*detDist * cosThetaVector[angI] + t * sinThetaVector[angI];

						//p2[0] = -1 * detDist * sinThetaVector[angI] + t * cosThetaVector[angI];   //Not needed for the optimized code
						//p2[1] = detDist * cosThetaVector[angI] + t * sinThetaVector[angI];

						//double ky = (p1[1]-p2[1])/(p1[0]-p2[0]);     //Optimized away using lookup table
						const double ky = minusCotanThetaVector[angI];      //This gives slightly different value for Theta=90deg which leads to small diff. compared to the non-optimized code
						//double pathInSinglePixel = sqrt(1+ky*ky)*pixSizes[0];   //Optimized away using lookup table
						const double pathInSinglePixel = pathInSinglePixelX[angI];

						//go through the columns of the image
						double yi_minus = ( halfPhantomHeight - (ky*( (colIdxMin-1)   *pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) * invPixSize1;
						const double yi_minusIncrement = -1*ky*invPixSize1*pixSizes[0];
						const double yi_plusIncrement = -1 * ky*pixSizeRatio01;
						for(int colIdx=colIdxMin; colIdx<=colIdxMax; ++colIdx){
							//double yi_minus = ( halfPhantomHeight - (ky*( colIdx   *pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) * invPixSize1;  //Always on the left side of the column
							//double yi_plus  = ( halfPhantomHeight - (ky*((colIdx+1)*pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) / pixSizes[1]; //Always on the right side of the column      //Optimized code in next line
							//Note: The heavy math wa pulled out from the for cycle. This leads slightly diferent results when Theta=90deg due to numerical errors
							yi_minus += yi_minusIncrement;
							double yi_plus = yi_minus + yi_plusIncrement;

							//Early exit criterion for performance
							if( (yi_minus < rowIdxMin-2) or (yi_plus < rowIdxMin-2) )
								continue;
							if( (yi_minus > rowIdxMax+2) or (yi_plus > rowIdxMax+2) )
								continue;

							//Early exit criterion for performance; This is more strict but the one above is faster
							/*if( (yi_minus <= rowIdxMin-1) and (yi_plus <= rowIdxMin-1) )
								continue;
							if( (yi_minus >= rowIdxMax+1) and (yi_plus >= rowIdxMax+1) )
								continue;*/

							int Yi_minusIdx, Yi_plusIdx;
							if( (yi_minus >= rowIdxMin) or (yi_plus >= rowIdxMin)){
								Yi_minusIdx = static_cast<int>(yi_minus+1)-1;     //This floor() is correct is yi_minus>-1 otherwise we don't care because it is out of bounds
								Yi_plusIdx = static_cast<int>(yi_plus+1)-1;
							}
							else{
								continue;
							}

							//double l_minus, l_plus; //intersecting lengths when crossing two pixels
							if( Yi_minusIdx == Yi_plusIdx ){ //intersecting only one pixel
								if( (Yi_minusIdx <= rowIdxMax) ){ // and (Yi_minusIdx >= rowIdxMin ) ){
									//l=sqrt(1+ky*ky)*pixSizes[0]; //Optimized away with pathInSinglePixel
									//sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += pathInSinglePixel * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
									sinogramPtr[sinogramPoint] += pathInSinglePixel * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
								}
							} else{
								if ( (Yi_minusIdx <= rowIdxMax) and (Yi_minusIdx >= rowIdxMin) ){
									const double l_minus=(std::max(Yi_minusIdx, Yi_plusIdx)-yi_minus) / (yi_plus - yi_minus) * pathInSinglePixel;
									sinogramPtr[sinogramPoint] += l_minus * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];

									//We have l_minus -> we can calculate l_plus with only a subtraction from pathInSinglePixel
									if( (Yi_plusIdx <= rowIdxMax) and (Yi_plusIdx >= rowIdxMin ) ){
										sinogramPtr[sinogramPoint] += (pathInSinglePixel- l_minus) * dataPtr[Yi_plusIdx*numberOfPixels[0] + colIdx];
										continue;
									}
								}

								if ( (Yi_plusIdx <= rowIdxMax) and (Yi_plusIdx >= rowIdxMin ) ){
									const double l_plus=(yi_plus - std::max(Yi_minusIdx, Yi_plusIdx)) / (yi_plus - yi_minus) * pathInSinglePixel;
									sinogramPtr[sinogramPoint] += l_plus * dataPtr[Yi_plusIdx*numberOfPixels[0] + colIdx];
								}
							}

						}
					}
				}
				else{      //beam intersects the rows at most two pixels
					//go through the different t values
					for(size_t pixIdx=static_cast<size_t>(minPixIdx); pixIdx<=static_cast<size_t>(maxPixIdx); ++pixIdx){
						const double t = pixPositions[pixIdx];

						size_t sinogramPoint = pixIdx + angI*static_cast<size_t>(pixNum);       //Index of the point in sinogram array

						p1[0]=detDist * sinThetaVector[angI] + t * cosThetaVector[angI];
						p1[1]=-1*detDist * cosThetaVector[angI] + t * sinThetaVector[angI];

						//p2[0] = -1 * detDist * sinThetaVector[angI] + t * cosThetaVector[angI]; //Not needed for the optimized code
						//p2[1] = detDist * cosThetaVector[angI] + t * sinThetaVector[angI];

						//double kx = (p1[0]-p2[0])/(p1[1]-p2[1]);   //Optimized away using lookup table
						double kx = minusTanThetaVector[angI];
						//double pathInSinglePixel = sqrt(1+kx*kx)*pixSizes[1];   //Optimized away using lookup table
						double pathInSinglePixel = pathInSinglePixelY[angI];

						//go through the rows of the image
						double xi_minus = (halfPhantomWidth + kx*( halfPhantomHeight - (rowIdxMin-1)     *pixSizes[1] - p1[1] ) + p1[0]  )  * invPixSize0;
						const double xi_minusIncrement = -1*kx*invPixSize0*pixSizes[1];
						const double xi_plusIncrement = -1* kx*pixSizeRatio10;
						for(int rowIdx=rowIdxMin; rowIdx<=rowIdxMax; ++rowIdx){
							//double xi_minusReal = (halfPhantomWidth + kx*( halfPhantomHeight - rowIdx     *pixSizes[1] - p1[1] ) + p1[0]  )  * invPixSize0;
							//double xi_plus  = (halfPhantomWidth + (kx*( halfPhantomHeight - (rowIdx+1) *pixSizes[1] - p1[1] ) + p1[0] ) ) / pixSizes[0];    //Optimized code in next line
							xi_minus += xi_minusIncrement;
							//double xi_plus = xi_minus - kx*pixSizeRatio10;
							double xi_plus = xi_minus + xi_plusIncrement;

							//Early exit criterion for performance
							if( (xi_minus < colIdxMin-2) or (xi_plus < colIdxMin-2) )
								continue;
							if( (xi_plus > colIdxMax+2) or (xi_plus > colIdxMax+2) )
								continue;

							//Early exit criterion for performance; This is more strict but the one above is faster
							/*if( (xi_minus <= colIdxMin-1) and (xi_plus <= colIdxMin-1) )
								continue;
							if( (xi_minus >= colIdxMax+1) and (xi_plus >= colIdxMax+1) )
								continue;*/

							int Xi_minusIdx, Xi_plusIdx;
							if( (xi_minus >= colIdxMin) or (xi_plus >= colIdxMin)){
								Xi_minusIdx = static_cast<int>(xi_minus+1)-1;
								Xi_plusIdx = static_cast<int>(xi_plus+1)-1;
							}
							else{
								continue;
							}

							double l_minus, l_plus; //intersecting lengths
							if( Xi_minusIdx == Xi_plusIdx ){
								if( (Xi_minusIdx <= colIdxMax) and (Xi_minusIdx >= colIdxMin ) ){
									//l=sqrt(1+kx*kx)*pixSizes[1]; //Optimized away with pathInSinglePixel
									sinogramPtr[sinogramPoint] += pathInSinglePixel * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
								}
							}
							else{
								if ( (Xi_minusIdx <= colIdxMax) and (Xi_minusIdx >= colIdxMin ) ){
									l_minus=(std::max(Xi_minusIdx, Xi_plusIdx)-xi_minus) / (xi_plus - xi_minus) * pathInSinglePixel;
									sinogramPtr[sinogramPoint] += l_minus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];

									//We have l_minus -> we can calculate l_plus with only a subtraction from pathInSinglePixel
									if ( (Xi_plusIdx <= colIdxMax) and (Xi_plusIdx >= colIdxMin ) ){   //Shortcut to avoid one more std::max call
										sinogramPtr[sinogramPoint] += (pathInSinglePixel-l_minus) * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
										continue;
									}
								}

								if ( (Xi_plusIdx <= colIdxMax) and (Xi_plusIdx >= colIdxMin ) ){
									l_plus=(xi_plus - std::max(Xi_minusIdx, Xi_plusIdx)) / (xi_plus - xi_minus) * pathInSinglePixel;
									sinogramPtr[sinogramPoint] += l_plus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
								}
							}
						}
					}
				}
			}
    	}
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Projection with ray-driven method ENDS here !!
    ////////////////////////////////////////////////////////////////////////////////

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedTime = stop - start;
    std::cout << "\nProjection with ray-driven method (optimized) took " << elapsedTime.count() << " milliseconds\n";
    return sinogram;
}

/***
 * Simple pixel-driven projector running on the CPU
 * @param actualPhantom The Phantom object which has to be projected
 * @param angles Projection angles [rad]
 * @return The sinogram matrix of the actualPhantom
 */
Eigen::MatrixXd Gen1CT::project_pixelDriven_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles){

	std::cout << "\nProjecting with pixel-driven method (interpolation) started";
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	auto pixSizes = actualPhantom.getPixSizes();
	auto numberOfPixels = actualPhantom.getNumberOfPixels();

	double t;
	const double piPer4 = M_PI/4;

	std::vector<double>	thetaVector,
	                    sinThetaVector, cosThetaVector,
	                    tanThetaVector, cotThetaVector,
	                    absSinThetaInvVector, absCosThetaInvVector;
	std::vector<bool> interpIsInY;

	for(int i=0; i<numAngles; i++){
		thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
		sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
		cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
		tanThetaVector.push_back( tan(thetaVector[static_cast<size_t>(i)]) );
		cotThetaVector.push_back( 1/tanThetaVector[static_cast<size_t>(i)] );
		absSinThetaInvVector.push_back( 1/std::abs(sinThetaVector[static_cast<size_t>(i)]) );
		absCosThetaInvVector.push_back( 1/std::abs(cosThetaVector[static_cast<size_t>(i)]) );

		interpIsInY.push_back( ( (thetaVector[static_cast<size_t>(i)] > piPer4 ) && (thetaVector[static_cast<size_t>(i)] < 3*piPer4) ) ||
				                 ( (thetaVector[static_cast<size_t>(i)] > 5*piPer4 ) && (thetaVector[static_cast<size_t>(i)] < 7*piPer4) ) );
	}

	for(size_t pixI=0; pixI<static_cast<size_t>(pixNum); ++pixI){
		t=pixPositions[pixI];
		for(size_t angI=0; angI < static_cast<size_t>(numAngles); ++angI){
			if( interpIsInY[angI] ){
				for(int objectXIndex=0; objectXIndex < numberOfPixels[0]; ++objectXIndex){
					double objectYinMM = t*sinThetaVector[angI]+ (t*cosThetaVector[angI] - actualPhantom.getXValueAtPix(objectXIndex))*cotThetaVector[angI];
					sinogram(static_cast<long>(pixI), static_cast<long>(angI) ) += actualPhantom.linear_atY(objectXIndex, objectYinMM) * absSinThetaInvVector[angI]*pixSizes[0];
				}
			} else{
				for(int objectYIndex=0; objectYIndex < numberOfPixels[1]; ++objectYIndex){
					double objectXinMM = t*cosThetaVector[angI] - (actualPhantom.getYValueAtPix(objectYIndex)-t*sinThetaVector[angI])*tanThetaVector[angI];
					sinogram(static_cast<long>(pixI), static_cast<long>(angI) ) += actualPhantom.linear_atX(objectYIndex, objectXinMM) * absCosThetaInvVector[angI]*pixSizes[1];
				}
			}
		}
	}

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "\nProjecting with pixel-driven method (interpolation) took " << elapsedTime.count() << " milliseconds" << std::endl;

	return sinogram;
}

/***
 * Projection with the improved Siddon method described in [1]
 * [1] Jacobs F et al. "A fast algorithm to calculate the exact radiological path through a pixel or voxel space
 * J. Comp. and Inf. Tech. - CIT 6, 1998, 1, p89-94
 * @param actualPhantom The Phantom object which has to be projected
 * @param angles Projection angles [rad]
 * @return The sinogram matrix of the actualPhantom
 */
Eigen::MatrixXd Gen1CT::project_Siddon_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles){

	std::cout << std::endl << "Projection with improved Siddon's algo started";
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();
	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	auto pixSizes = actualPhantom.getPixSizes();
	auto numberOfPixels = actualPhantom.getNumberOfPixels();

	std::vector<double>	thetaVector,
						sinThetaVector, cosThetaVector;

	for(int i=0; i<numAngles; i++){
		thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
		sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
		cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
	}

	const double* dataPtr=actualPhantom.getDataAsEigenMatrixRef().data();

	double detDist = 1.1 * sqrt(pow(pixSizes[0]*numberOfPixels[0], 2) + pow(pixSizes[1]*numberOfPixels[1], 2) ); //Distance between the detector plane and centre of rotation

	double dX =actualPhantom.getPixSizes()[0];
	double dY =-1*actualPhantom.getPixSizes()[1];

	double bX = actualPhantom.getXValueAtPix(0) - dX/2;
	double bY = actualPhantom.getYValueAtPix(0) - dY/2;

	int Nx = actualPhantom.getNumberOfPixels()[0]+1;
	int Ny = actualPhantom.getNumberOfPixels()[1]+1;

	int iMin, iMax, jMin, jMax;

	//Improved Siddon's Algo:
	for(size_t angI=0; angI < static_cast<size_t>(numAngles); ++angI){
		for(size_t pixI=0; pixI < static_cast<size_t>(pixNum); ++pixI){
			double t=pixPositions[pixI];
			std::array<double,2> p1{ detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
									 -1*detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};
			std::array<double,2> p2{ -1 * detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
									 detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};

			double xDiff = p2[0] - p1[0];
			double yDiff = p2[1] - p1[1];

			double d_conv = sqrt(xDiff*xDiff + yDiff*yDiff);

			double d12 = 0; //Accumulator for the raySum

			if((angles[static_cast<long>(angI)] == 0) || angles[static_cast<long>(angI)] == M_PI){ //Edge case when ray is vertical
				int colIndex=std::floor( (p1[0]-bX)/dX );
				if (colIndex >=0 && colIndex < Nx-1){
					//Sum the matrix in vertical direction
					for(int rowIndex=0; rowIndex < Ny-1; ++rowIndex){
						d12 += dataPtr[ rowIndex*(Nx-1)+ colIndex];
					}
				}
				d12 *= -1*dY;
			}
			else if ((angles[static_cast<long>(angI)] == M_PI/2) || angles[static_cast<long>(angI)] == 3*M_PI/2){ //Edge case when ray is horizontal
				int rowIndex=std::floor( (p1[1]-bY)/dY );
				if (rowIndex >=0 && rowIndex < Ny-1){
					//Sum the matrix in vertical direction
					for(int colIndex=0; colIndex < Nx-1; ++colIndex){
						d12 += dataPtr[ rowIndex*(Nx-1)+ colIndex];
					}
				}
				d12 *= dX;
			}
			else{  //General case

				double alp_xu = dX / std::abs(p2[0] - p1[0]); //Update of alpha when step in X
				double alp_yu = std::abs(dY) / std::abs(p2[1] - p1[1]); //Update of alpha when step in Y

				double i_u = (p1[0]<p2[0]) ? 1 : -1;  // Update of index i when stepping in X
				double j_u = (p1[1]<p2[1]) ? -1 : 1;  // Update of index j when stepping in Y

				double alp_min; //Entry point
				double alp_max; //Exit point

				double alp_xmin, alp_xmax;
				double alp_ymin, alp_ymax;

				double alp_xBegin, alp_xEnd; //Refers to alpha_x(0) and alpha_x(N_x-1) in the article
				double alp_yBegin, alp_yEnd; //Refers to alpha_y(0) and alpha_y(N_x-1) in the article

				alp_xBegin = (bX-p1[0])/xDiff;  // \alpha when entering 0th index of pixel space in X direction
				alp_yBegin = (bY-p1[1])/yDiff;  // \alpha when entering 0thindex of pixel space in Y direction

				alp_xEnd = (bX+(Nx-1)*dX-p1[0])/xDiff;  // \alpha when entering last index of pixel space in X direction
				alp_yEnd = (bY+(Ny-1)*dY-p1[1])/yDiff;  // \alpha when entering last index of pixel space in Y direction

				alp_xmin = std::min(alp_xBegin, alp_xEnd);
				alp_xmax = std::max(alp_xBegin, alp_xEnd);
				alp_ymin = std::min(alp_yBegin, alp_yEnd);
				alp_ymax = std::max(alp_yBegin, alp_yEnd);

				alp_min = std::max(alp_xmin, alp_ymin); // \alpha when entering pixel space
				alp_max = std::min(alp_xmax, alp_ymax); // \alpha when entering pixel space

				if (alp_min > alp_max){   //The ray does not intersect the pixel space
					continue;
				}

				if(p1[0] <p2[0]){
					//going from left to right
					if( alp_min == alp_xmin){  // Entering pixel space through the left side
						iMin = 1;
					} else{
						iMin = std::ceil( (p1[0]+alp_min*(p2[0]-p1[0]) - bX) / dX );  //Entering in a general position
					}
					if( alp_max == alp_xmax){ // Leaving pixel space on the right side
						iMax = Nx-1;
					} else {
						iMax = std::floor( (p1[0]+alp_max*(p2[0]-p1[0]) - bX) / dX ); //Leaving in a general position
					}
				}
				else{
					//going from right to left
					if( alp_min == alp_xmin){
						iMax = Nx-2;   // Entering pixel space through the right side
					} else{
						iMax = std::floor( (p1[0]+alp_min*(p2[0]-p1[0]) - bX) / dX ); //Entering in a general position
					}
					if( alp_max == alp_xmax){ // Leaving pixel space on the left side
						iMin = 0;
					} else {                  //Leaving in a general position
						iMin = std::ceil ( (p1[0]+alp_max*(p2[0]-p1[0]) - bX) / dX );
					}
				}

				if(p1[1] < p2[1]){
					//going from bottom to top
					if( alp_min == alp_ymin){  // Entering pixel space through the bottom side
						jMax = Ny - 2;
					} else{
						jMax = std::floor( (p1[1]+alp_min*(p2[1]-p1[1]) - bY) / dY );  //Entering in a general position
					}
					if( alp_max == alp_ymax){ // Leaving pixel space through the top
						jMin = 0;
					} else {
						jMin = std::ceil( (p1[1]+alp_max*(p2[1]-p1[1]) - bY) / dY );  //Leaving in a general position
					}
				}
				else{
					//going from top to bottom
					if( alp_min == alp_ymin){ // Entering pixel space through the top side
						jMin = 1;
					} else{
						jMin = std::ceil( (p1[1]+alp_min*(p2[1]-p1[1]) - bY) / dY ); //Entering in a general position
					}
					if( alp_max == alp_ymax){ // Leaving pixel space through the bottom
						jMax = Ny-1;
					} else {
						jMax = std::floor( (p1[1]+alp_max*(p2[1]-p1[1]) - bY) / dY );//Leaving in a general position
					}
				}

				double alpX;
				if (p1[0] < p2[0]) //We are going from left to right
					alpX = ( bX + iMin*dX - p1[0]) / xDiff;
				else               //We are going from right to left
					alpX = ( bX + iMax*dX - p1[0]) / xDiff;

				double alpY;
				if (p1[1] < p2[1]) //We are going from bottom to top
					alpY = ( bY + jMax*dY - p1[1]) / yDiff;
				else               //We are going from top to bottom
					alpY = ( bY + jMin*dY - p1[1]) / yDiff;

				int Np = (iMax-iMin+1) + (jMax-jMin+1);

				double argument = (std::min(alpX, alpY) + alp_min)/2;
				int i = std::floor( (p1[0] + argument*(p2[0]-p1[0]) - bX) / dX );  //Index of the first intersected pixel
				int j = std::floor( (p1[1] + argument*(p2[1]-p1[1]) - bY) / dY );  //Index of the first intersected pixel

				double alp_c= alp_min;

				for(int counter =0; counter<Np; counter++){
					if(alpX < alpY){
						d12 += (alpX-alp_c) * dataPtr[ j*(Nx-1) + i];
						i += i_u;
						alp_c=alpX;
						alpX += alp_xu;
					} else{
						d12 += (alpY-alp_c) * dataPtr[ j*(Nx-1) + i];
						j += j_u;
						alp_c=alpY;
						alpY += alp_yu;
					}
				}
				d12 *= d_conv;
			}

			//Update the sinogram
			sinogram(static_cast<long>(pixI), static_cast<long>(angI)) += d12;
		}
	}

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "\nProjection with Siddon's algorithm took " << elapsedTime.count() << " milliseconds" << std::endl;

	return sinogram;
}

/***
 * Get a copy of a measured sinogram from the stored library
 * @param label Label string of the sinogram
 * @return Copy of the sinogram if the label found, othervise an empty CTScan object with "EMPTY_SCAN" label
 */
CTScan Gen1CT::getMeasurement(const std::string& label){
	if(scans.find(label) != scans.end()){
			return scans.at(label);
	}
	else{
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! returning an empty CTScan!!";
		return CTScan("EMPTY_SCAN", Eigen::MatrixXd::Zero(pixNum, 1), detWidth, Eigen::VectorXd{0}, I0);
	}
}

/***
 * Display the required sinogram from the library in a separate window
 * @param label  Label string of the sinogram
 */
void Gen1CT::displayMeasurement(const std::string& label){
	if(scans.find(label) != scans.end()){
			scans.at(label).display(label);
	}
	else
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

/***
 * Backproject a CTScan object using a selected backprojector
 * @param sinogram Sinogram which needs to be backprojected
 * @param numberOfRecPoints Number or backprojection points along the X and Y axis
 * @param resolution Distance between the reconstruction points
 * @param backProjector Type of backprojector
 * @return Matrix with the backprojected values
 */
Eigen::MatrixXd Gen1CT::backProject(const CTScan& sinogram,
									const std::array<int,2>& numberOfRecPoints,
									const std::array<double,2>& resolution,
									backprojectorType backProjector){
	Eigen::MatrixXd backProjection;

	switch(backProjector){
		case backprojectorType::pixelDriven:
			backProjection = backProject_pixelDriven_CPU(sinogram, numberOfRecPoints, resolution);
			break;
		case backprojectorType::rayDriven:
			backProjection = backProject_rayDriven_CPU(sinogram, numberOfRecPoints, resolution);
			break;
#if ENABLE_CUDA
		case backprojectorType::rayDriven_GPU:
			backProjection = backProject_rayDriven_GPU(sinogram, numberOfRecPoints, resolution);
			break;
#endif
	}
	return backProjection;
}

/**
 * Backproject interpolated values of the sinogram to the pixels
 *
 * @param sinogram CTScan object which has to be backprojected
 * @param numberOfRecPoints Number of points in the backprojection
 * @param resolution Pixel size of the backprojected image
 * @return Matrix with the backprojected values
 */
Eigen::MatrixXd Gen1CT::backProject_pixelDriven_CPU(const CTScan& sinogram,
						 const std::array<int,2>& numberOfRecPoints,
		                 const std::array<double,2>& resolution){

	std::cout << "Backprojection with pixel-driven method started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	assert(numberOfRecPoints.size()==2);
	assert(resolution.size()==2);

	Eigen::MatrixXd backprojection = Eigen::MatrixXd::Zero(numberOfRecPoints[0], numberOfRecPoints[1]);

	//Vectors with coordinates of the grid in real space
	double xMax=numberOfRecPoints[0]*resolution[0]/2;
	double xMin=-1*xMax;
	Eigen::VectorXd xValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[0], xMin+resolution[0]/2, xMax-resolution[0]/2)};

	double yMax=numberOfRecPoints[1]*resolution[1]/2;
	double yMin=-1*yMax;
	Eigen::VectorXd yValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[1], yMax-resolution[1]/2, yMin+resolution[1]/2)};

	const Eigen::VectorXd& angs=sinogram.getAnglesConstRef();
	const Eigen::MatrixXd& sinoData=sinogram.getDataAsEigenMatrixRef();

	std::vector<double> sinTheta(static_cast<size_t>(angs.size()));
	std::vector<double> cosTheta(static_cast<size_t>(angs.size()));
	for(unsigned int i=0; i<angs.size(); ++i){
		sinTheta[i]=sin(angs[i]);
		cosTheta[i]=cos(angs[i]);
	}

	//For each point in real space
	double offset = detWidth/2 - detWidth/pixNum/2;
	double invPixRes = pixNum/detWidth;
	double minPixPosition = -1*offset;
	double maxPixPosition = offset;
	//For every pixel
	for(int xIdx=0; xIdx<numberOfRecPoints[0]; ++xIdx){
		for(int yIdx=0; yIdx<numberOfRecPoints[1]; ++yIdx){
			//For every angle
			for(unsigned int thIdx=0; thIdx<angs.size(); ++thIdx){
				//Add the corresponding interpolated points from the sinogram
				double tValue = xValues[xIdx]*cosTheta[thIdx] + yValues[yIdx]*sinTheta[thIdx];
				if( (tValue<minPixPosition) || (tValue > maxPixPosition))
					continue;
				double pixIdxDouble = (tValue + offset) * invPixRes;  //pixel coordinate
				int pixIdxInt = round(pixIdxDouble);     //pixel bin index

				double valueInLowerPixel, valueInHigherPixel;
				if(pixIdxInt > pixIdxDouble ){
					valueInLowerPixel  = sinoData(pixIdxInt-1, thIdx);
					valueInHigherPixel = sinoData(pixIdxInt, thIdx);
					backprojection(xIdx, yIdx) += valueInLowerPixel + (valueInHigherPixel - valueInLowerPixel) *
							                                           (1-(pixIdxInt-pixIdxDouble));
				}
				else{
					valueInLowerPixel  = sinoData(pixIdxInt, thIdx);
					valueInHigherPixel = sinoData(pixIdxInt+1, thIdx);
					backprojection(xIdx, yIdx) += valueInLowerPixel + (valueInHigherPixel - valueInLowerPixel) *
							                                           (pixIdxDouble-pixIdxInt);
				}
			}
		}
	}
	//Multiply with dTheta
	backprojection = backprojection*M_PI/angs.size();

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "Backprojection with pixel-driven method took " << elapsedTime.count() << " milliseconds";

	return backprojection;
}

/**
 * Backproject using the method proposed by Hao Gao in [1]
 * [1] Med. Phys. 39(11), 7110 (2012)
 *
 * @param sinogram  The sinogram which is backprojected
 * @param numberOfRecPoints Number of reconstruction points
 * @param resolution Resolution in [mm]
 * @return Eigen::MatrixXd with the reconstructed data
 */
Eigen::MatrixXd Gen1CT::backProject_rayDriven_CPU(const CTScan& sinogram,
						 const std::array<int,2>& numberOfRecPoints,
		                 const std::array<double,2>& resolution){

	std::cout << "Backprojection using ray-driven method started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	assert(numberOfRecPoints.size()==2);
	assert(resolution.size()==2);

	Eigen::MatrixXd backProjection = Eigen::MatrixXd::Zero(numberOfRecPoints[0], numberOfRecPoints[1]);

	//Vectors with coordinates of the grid in real space
	double xMax=numberOfRecPoints[0]*resolution[0]/2;
	double xMin=-1*xMax;
	Eigen::VectorXd xValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[0], xMin+resolution[0]/2, xMax-resolution[0]/2)};

	double yMax=numberOfRecPoints[1]*resolution[1]/2;
	double yMin=-1*yMax;
	Eigen::VectorXd yValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[1], yMax-resolution[1]/2, yMin+resolution[1]/2)};


	const Eigen::VectorXd& angles=sinogram.getAnglesConstRef();
	int numAngles = angles.size();
	const Eigen::MatrixXd& sinoData=sinogram.getDataAsEigenMatrixRef();

	////////////////////////////////////////////////////////
	/// START the method of Hao Gao
	////////////////////////////////////////////////////////

    std::vector<double>	thetaVector,
    	                sinThetaVector,
						cosThetaVector;

    for(int i=0; i<numAngles; i++){
    	thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
    	sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
    	cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
    }

    const double pixelRadius = sqrt(pow(resolution[0]/2,2) + pow(resolution[1]/2,2)); //Radius of circle drawn around a pixel
    const double halfDetWidth = detWidth / 2;
    const double invDetPixSize= pixNum/detWidth;
	//Go through the pixels
	for(int xIdx=0; xIdx<numberOfRecPoints[0]; ++xIdx){
		for(int yIdx=0; yIdx<numberOfRecPoints[1]; ++yIdx){
			//Go through the angles
			for(int angIdx=0; angIdx<numAngles; ++angIdx){
				//Determine the contributing detector pixels (i.e. rays)
				double xr = xValues(static_cast<long>(xIdx))*cosThetaVector[static_cast<size_t>(angIdx)]+
						       yValues(static_cast<long>(yIdx))*sinThetaVector[static_cast<size_t>(angIdx)];   //This corresponds to "t" of the image pixel center

				double minPixIdx = (xr + halfDetWidth - pixelRadius )*invDetPixSize;
				double maxPixIdx = (xr + halfDetWidth + pixelRadius )*invDetPixSize;

				//Calculate the intersection length and add to the backprojected image
				//Go through the possible pixels
				double lSum=0;
				double angleContrib = 0;
				for(int detPixIdx=std::max(0,static_cast<int>(minPixIdx));
					    detPixIdx <= std::min(static_cast<int>(maxPixIdx)+0, static_cast<int>(pixNum-1) );
					    ++detPixIdx){
					double a = cosThetaVector[static_cast<size_t>(angIdx)];
					double b = sinThetaVector[static_cast<size_t>(angIdx)];
					double c = pixPositions[static_cast<size_t>(detPixIdx)];    // "t" of the detector pixel center

					//Calculate the boundaries of pixel intersection types
					double d_max = (std::abs(a*resolution[0]) + std::abs(b*resolution[1]))/2; //For definition see HaoGao's article
					double d_min =  std::abs( std::abs(a*resolution[0]) - std::abs(b*resolution[1]) )/2;

					//Calculate the intersection length
					double d_act=std::abs(a*xValues(xIdx) + b*yValues(yIdx) - c);
					double l;

					if(d_act < d_min){
						if(std::abs(a) < std::abs(b)){
							l=resolution[0]/std::abs(b);
						}
						else{
							l=resolution[1]/std::abs(a);
						}
					}
					else if( (d_act >= d_min) and (d_act < d_max) ){
						l=(d_max-d_act)/std::abs(a)/std::abs(b);
					}
					else{
						l=0;
					}

					lSum += l;
					angleContrib += sinoData(detPixIdx, angIdx)*l;
				}
				if(lSum != 0){
					backProjection(xIdx, yIdx) += angleContrib/lSum;
				}
			}
		}
	}
	//Multiply with dTheta
	backProjection = backProjection*M_PI/numAngles;

	////////////////////////////////////////////////////////
	/// END of the method of Hao Gao
	////////////////////////////////////////////////////////

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "Backprojection using ray-driven method took " << elapsedTime.count() << " milliseconds" << std::endl;

	return backProjection;
}

/***
 * Filtering using the Filter functor
 * @param sinogramID Label of the sinogram to be filtered
 * @param filter Filter class object
 * @return CTScan object with the filtered sinogram or an empty CTScan if sinogramID was not found
 */
CTScan Gen1CT::applyFilter(const std::string& sinogramID, Filter filter){

	//Timing
	std::cout << "Filtering started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	//Getting the reference to the sinogram to be filtered
	if(scans.find(sinogramID) == scans.end()){
				std::cout << "\n Sinogram was not found, exiting from filtering function";
				return CTScan();
	}
	CTScan& sinogram = scans.at(sinogramID);
	int numAngles = sinogram.getAnglesConstRef().rows();

	//Zero padding of the sinogram
	int pixNumPadded = std::pow(2, std::ceil(std::log2(2*pixNum-1)));
	std::cout << "\nOriginal pixNum:" << pixNum << " ZeroPadded to: " << pixNumPadded <<"\n";

	Eigen::MatrixXd paddedSinogram = Eigen::MatrixXd::Zero(pixNumPadded, numAngles );
	int startIndex=floor((pixNumPadded-pixNum)/2);
	paddedSinogram.block(startIndex, 0, pixNum, numAngles) = sinogram.getDataAsEigenMatrixRef();

	//Fourier transform of the padded Sinogram:
	Eigen::FFT<double> fft;
	Eigen::MatrixXcd fftOfSinogram = Eigen::MatrixXcd::Zero(pixNumPadded, numAngles);
	for(int i=0; i<paddedSinogram.cols(); i++){
		fftOfSinogram.col(i) = fft.fwd(paddedSinogram.col(i));
	}

	filter(fftOfSinogram);

	//IFFT of filtered sinogram
	double tau=detWidth/pixNum;
	for(int i=0; i<fftOfSinogram.cols(); i++){
		paddedSinogram.col(i) = 1/tau * fft.inv(fftOfSinogram.col(i)).real();
	}

	auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "Filtering took " << elapsedTime.count() << " milliseconds" << std::endl;

	return CTScan(sinogramID+"F", paddedSinogram.block(startIndex,0, pixNum, numAngles),
		detWidth, sinogram.getAnglesConstRef(), sinogram.getI0()) ;
}

/***
 * Reconstruct an image with Filtered BackProjection method
 * @param sinogramID The sinogram which has to be reconstructed
 * @param numberOfRecPoints Number of reconstruction points in X and Y directions
 * @param resolution Sepsize in X and Y directions
 * @param filterType Type of the filter used for reconstruction
 * @param cutOffFreq Cut off frequency [0-1]
 * @param backProjectAlgo Type of backprojector
 * @param imageID Label of the reconstructed image which will be used to store data
 * @param referenceImage Reference image for calculation of L2 norm of error (optional)
 */
void Gen1CT::filteredBackProject(std::string sinogramID,
			                    const std::array<int,2>& numberOfRecPoints,
								const std::array<double,2>& resolution,
								FilterType filterType,
								double cutOffFreq,
								backprojectorType backProjectAlgo,
								const std::string& imageID,
								std::string referenceImage){
	//Timing
    std::cout << "\nFiltered Backprojection started";
	auto start = std::chrono::high_resolution_clock::now();

	if(scans.find(sinogramID) == scans.end()){
			std::cout << std::endl << "ERROR!! sinogramID: \"" << sinogramID << "\" could not be found!! Abort mission";
			return;
	}

	//Convert the counts to line integrals
	CTScan tmpFBPSinogram = scans[sinogramID];
	tmpFBPSinogram.convertToLineIntegrals();
	scans.emplace("tmpSinogram", tmpFBPSinogram);
	sinogramID="tmpSinogram";

	//Apply filter on the sinogram
	CTScan filteredScan = applyFilter(sinogramID, Filter(filterType, cutOffFreq) );

	//Backproject the image
	Reconst recImage(imageID, backProject(filteredScan, numberOfRecPoints,resolution, backProjectAlgo), resolution);

	//Calculate the L2 norm of difference to reference image
	double normError=0.0;
	if(referenceImage != ""){
		if(phantoms.find(referenceImage) == phantoms.end()){
			std::cout << std::endl << "ERROR!! referenceImage: \"" << referenceImage << "\" could not be found!! L2norm could not be calculated!";
		} else{
			normError = recImage.compareNorm(phantoms[referenceImage]);
			std::cout << "\nL2 normalized error: " << normError;
		}
	}

	//Move the backprojected image to reconsts map
	auto it = reconsts.find(imageID);
	if(it != reconsts.end()){
		std::cout << std::endl << "WARNING! A reconstructed image with label \"" << imageID << "\" already exists!!! Overwriting!!!";
		reconsts.erase(it);
	}

	auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "\nReconstruction with Filtered BackProjection took " << elapsedTime.count()/1000.0 << " seconds" << std::endl;

	reconsts.emplace(imageID, Reconst(imageID, recImage.getDataAsEigenMatrixRef(), resolution, std::vector<double>{normError}));
	scans.erase("tmpSinogram");
}

/***
 * Display a refonstruction from the library of reconstructions
 * @param label The label of the reconstruction to be shown
 */
void Gen1CT::displayReconstruction(const std::string& label){
	if(reconsts.find(label) != reconsts.end()){
		reconsts[label].display(label);
	}else
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

/***
 * Run reconstruction with the emission ML-EM method. This method is not ideal for low-count measurements.
 * @param sinogramID Label of the sinogram to be backprojected
 * @param numberOfRecPoints Number of reconstruction points in X and Y direction
 * @param resolution Distance between reconstruction points in X and Y direction
 * @param projectAlgo Type of projector used for the projection step
 * @param backprojectAlgo Type of projector used for the backprojection step
 * @param imageID Label for the reconstructed image
 * @param numberOfIterations Number of iteration steps performed
 * @param referenceImage Reference image for error calculation
 */
void Gen1CT::MLEMReconst(std::string sinogramID,
			             const std::array<int,2>& numberOfRecPoints,
					     const std::array<double,2>& resolution,
						 projectorType projectAlgo,
						 backprojectorType backprojectAlgo,
						 const std::string& imageID,
						 int numberOfIterations,
						 std::string referenceImage){
	//Timing
	std::cout << "Reconstruction with emission ML-EM started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	if(scans.find(sinogramID) == scans.end()){
				std::cout << std::endl << "ERROR!! sinogramID: \"" << sinogramID << "\" could not be found!! Abort mission";
				return;
	}
	CTScan actualScan = scans.at(sinogramID);

	//Convert from counts to line integrals.
	//This is only an APPROXIMATION!! -> Not valid for low counts!!!
	actualScan.convertToLineIntegrals();

	//Sinorgam consisting of const ones for the normalization factor calculation
	CTScan constOnes("constOnes",
			         Eigen::MatrixXd::Ones(pixNum, actualScan.getNumberOfPixels()[1] ),
					 detWidth,
					 actualScan.getAnglesConstRef(),
					 0.0); //I0 is Zero

	//Calculate the normalization factors
	Phantom normFactors;
	normFactors = Phantom("normFactors",
		      	          backProject(constOnes, numberOfRecPoints, resolution, backprojectAlgo),
						  resolution);

	//Initialize the reconstructed image with const ones
	Phantom reconstImage("reconstructedImage",
			             Eigen::MatrixXd::Ones(numberOfRecPoints[0], numberOfRecPoints[1])*muWater,
						 resolution);

	std::vector<double> differenceNorms(0);
	for(int itNumber = 0; itNumber < numberOfIterations; ++itNumber){

		std::cout << "\nIteration:" << itNumber+1 << " / " << numberOfIterations;

		//Forward project the estimation
		CTScan forwardProj( "iteration",
							project(reconstImage, actualScan.getAnglesConstRef(), projectAlgo),
							detWidth,
					        actualScan.getAnglesConstRef(),
							0.0); //I0=0 because the Poisson statistics is not applied

		//Create the correction image
		//backproject the measurement elementwise divided with forwardprojection
		Phantom corrImage;
		corrImage = Phantom("corrImage",
				           backProject( actualScan/(forwardProj+0.001),
								                   numberOfRecPoints,
				                                   resolution,
												   backprojectAlgo),
		                   resolution);

		reconstImage = reconstImage / normFactors * corrImage;

		if(referenceImage != ""){
			if(phantoms.find(referenceImage) == phantoms.end()){
				std::cout << std::endl << "ERROR!! referenceImage: \"" << referenceImage << "\" could not be found!! L2norm could not be calculated!";
				differenceNorms.push_back(0.0);
			} else{
				differenceNorms.push_back(reconstImage.compareNorm(phantoms[referenceImage]));
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "Reconstruction with transmission ML-EM took " << elapsedTime.count()/1000.0 << " seconds" << std::endl;

	//Move the backprojected image to reconsts map
	auto it = reconsts.find(imageID);
	if(it != reconsts.end()){
		std::cout << std::endl << "WARNING! A reconstructed image with label \"" << imageID << "\" already exists!!! Overwriting!!!";
		reconsts.erase(it);
	}
	reconsts.emplace(imageID, Reconst(imageID, reconstImage.getDataAsEigenMatrixRef(), resolution, differenceNorms));
}

/***
 * Compare the same horizontal or vertical cuts (rows or columns)
 * of the phantom and the reconstruction
 *
 * @param direction  'Y' or 'X' for rows or columns respectively
 * @param position   position in [mm] of the comparision
 * @param phantomID  ID of the phantom
 * @param reconstID  ID of the reconstruction
 */
void Gen1CT::compareRowPhantomAndReconst(char direction, double position,
		                                 const std::string& phantomID, const std::string& reconstID){

	if(phantoms.find(phantomID) == phantoms.end()){
		std::cout << std::endl << "ERROR!! phantomID: \"" << phantomID << "\" could not be found!! Aborting compareRowPhantomAndReconst function";
		return;
	}

	if(reconsts.find(reconstID) == reconsts.end()){
		std::cout << std::endl << "ERROR!! reconstID: \"" << reconstID << "\" could not be found!! Aborting compareRowPhantomAndReconst function";
		return;
	}

	std::array<double,2> phantomPixSizes= phantoms[phantomID].getPixSizes();
	std::array<int,2> phantomPixNum  = phantoms[phantomID].getNumberOfPixels();

	std::array<double,2> recPixSizes=reconsts[reconstID].getPixSizes();
	std::array<int,2> recPixNum  =reconsts[reconstID].getNumberOfPixels();

	if( direction == 'Y' ){ //cut at the y=position line
		int recRowNum = std::round( (recPixSizes[1]*recPixNum[1]/2 - position)/recPixSizes[1] );
		if ( (recRowNum < 0) or (recRowNum >= recPixNum[1]) ){
			std::cout << "Not possible to show the required slice of the reconstructed image!!";
			return;
		}
		int phantomRowNum = std::round( (phantomPixSizes[1]*phantomPixNum[1]/2 - position)/phantomPixSizes[1] );
		if ( (phantomRowNum < 0) or (phantomRowNum >= phantomPixNum[1]) ){
			std::cout << "Not possible to show the required slice of the phantom!!";
			return;
		}

		Eigen::VectorXd BPSlice = reconsts[reconstID].getDataAsEigenMatrixRef().col(recRowNum);
		Eigen::VectorXd ObjSlice = phantoms[phantomID].getDataAsEigenMatrixRef().col(phantomRowNum);

		std::vector<float> recXvals(static_cast<size_t>(recPixNum[1]));
		for(int idx=0; idx<recPixNum[1]; ++idx){
			recXvals[static_cast<size_t>(idx)] = -1*recPixSizes[1]*recPixNum[1]/2 + (idx+0.5)*recPixSizes[1];
		}
		std::vector<float> phantomXvals(static_cast<size_t>(phantomPixNum[1]));
		for(int idx=0; idx<phantomPixNum[1]; ++idx){
			phantomXvals[static_cast<size_t>(idx)] = -1*phantomPixSizes[1]*phantomPixNum[1]/2 + (idx+0.5)*phantomPixSizes[1];
		}

		auto h=matplot::figure();
		matplot::plot(recXvals, std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
		matplot::hold(true);
		matplot::plot(phantomXvals, std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
		h->show();
	}
	else if (direction == 'X'){
		int recColNum = std::round( (position + recPixSizes[0]*recPixNum[0]/2)/recPixSizes[0] );
		if ( (recColNum < 0) or (recColNum >= recPixNum[0]) ){
			std::cout << "Not possible to show the required slice of the reconstructed image!!";
			return;
		}
		int phantomColNum = std::round( (position + phantomPixSizes[0]*phantomPixNum[0]/2)/phantomPixSizes[0] );
		if ( (phantomColNum < 0) or (phantomColNum >= phantomPixNum[0]) ){
			std::cout << "Not possible to show the required slice of the phantom!!";
			return;
		}
		Eigen::VectorXd BPSlice = reconsts[reconstID].getDataAsEigenMatrixRef().row(recColNum);
		Eigen::VectorXd ObjSlice = phantoms[phantomID].getDataAsEigenMatrixRef().row(phantomColNum);

		std::vector<float> recYvals(static_cast<size_t>(recPixNum[0]));
		for(int idx=0; idx<recPixNum[0]; ++idx){
			recYvals[static_cast<size_t>(idx)] = recPixSizes[0]*recPixNum[0]/2 - (idx+0.5)*recPixSizes[0];
		}
		std::vector<float> phantomYvals(static_cast<size_t>(phantomPixNum[0]));
		for(int idx=0; idx<phantomPixNum[0]; ++idx){
			phantomYvals[static_cast<size_t>(idx)] = phantomPixSizes[0]*phantomPixNum[0]/2 - (idx+0.5)*phantomPixSizes[0];
		}

		auto h=matplot::figure();
		matplot::plot(recYvals, std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
		matplot::hold(true);
		matplot::plot(phantomYvals, std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
		h->show();
	}
	else{
		std::cout << "Unable to determine the slicing direction in compareRowPhantomAndReconst function."
				  << "\n direction should be \'X\' or \'Y\'. Aborting function";
		return;
	}
}

/***
 * Print the parameters of a Phantom from the library
 * @param phantomLabel Label of the phantom of which the parameters has to be printed
 */
void Gen1CT::printPhantomParams(const std::string& phantomLabel){
	std::cout << "\n Phantom name: " << phantomLabel;
	std::cout << "\n Number of pixels: " << phantoms[phantomLabel].getNumberOfPixels()[0] << " x " << phantoms[phantomLabel].getNumberOfPixels()[1];
	std::cout << "\n Pixel sizes: " << phantoms[phantomLabel].getPixSizes()[0] << " x " << phantoms[phantomLabel].getPixSizes()[1];
}

/***
 * Reconstruct the sinogram using the Separable Paraboloidal Surrogates (SPS) algorithm
 * References: [1] H. Erdogan and J. A. Fessler. Ordered subsets algorithms for transmission tomography.
 *                    Phys. Med. Biol., 44(11):2835–51, November 1999.
 *             [2] J. A. Fessler. Grouped coordinate descent algorithms for robust edge-preserving image restoration.
 *                    Proc. SPIE 3170 Im. Recon. and Restor. II, pages 184–94, 1997.
 *             [3] J. Fessler Statistical Image Reconstruction MEthods for Transmission Tomography
 *                    SPIE Handbook of Medical Imaging, Vol I.pp 1-70, 2000
 *
 * @param sinogramID Label of the sinogram to be backprojected
 * @param numberOfRecPoints Number of reconstruction points in X and Y direction
 * @param resolution Distance between reconstruction points in X and Y direction
 * @param projectAlgo Type of projector used for the projection step
 * @param backProjectAlgo Type of projector used for the backprojection step
 * @param imageID Label for the reconstructed image
 * @param numberOfIterations Number of iteration steps performed
 * @param regularizerFunction Type of the regularization function
 * @param beta Regularization parameter
 * @param delta Regularization parameter
 * @param referenceImage Reference image for error calculation
 */
void Gen1CT::SPSReconst(std::string sinogramID,
						const std::array<int,2>& numberOfRecPoints,
						const std::array<double,2>& resolution,
						projectorType projectAlgo,
						backprojectorType backProjectAlgo,
						const std::string& imageID,
						int numberOfIterations,
						regularizerType regularizerFunction,
						double beta,
						double delta,
						std::string referenceImage){
	//Timing
    std::cout << "Reconstruction with SPS method started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	if(scans.find(sinogramID) == scans.end()){
		std::cout << std::endl << "ERROR!! sinogramID: \"" << sinogramID << "\" could not be found!! Abort mission";
		return;
	}
	CTScan actualScan = scans.at(sinogramID);

	//Phantom of const numbers:
	Phantom constOnesPhantom("constOnesPhantom",
			         Eigen::MatrixXd::Ones(numberOfRecPoints[0], numberOfRecPoints[1] ),
					 resolution);

	//Calculate the normalization factors
	CTScan normFactors("normFactors",
		      	       project(constOnesPhantom, actualScan.getAnglesConstRef(), projectAlgo),
					   detWidth,
					   actualScan.getAnglesConstRef(),
					   0.0); //I0=0 because the Poisson statistics is not applied

	//Iteration converge slowly -> start with the FBP image
	filteredBackProject(sinogramID, numberOfRecPoints,
			resolution, FilterType::RamLak, 0.5, backProjectAlgo,
			"FBPrecImage");

	Phantom reconstImage("reconstructedImage",
			             reconsts["FBPrecImage"].getDataAsEigenMatrixRef().cwiseMax(0.0),
						 resolution);
	reconsts.erase("FBPrecImage");

	std::vector<double> differenceNorms;
	//START the ITERATION
	for(int itNumber = 0; itNumber < numberOfIterations; ++itNumber){

		std::cout << "\nIteration:" << itNumber+1 << " / " << numberOfIterations;

		//Calculate l_i^(n) = [A\mu^(n)]_i Forward-project the actual estimate
		CTScan li("li",
				  project(reconstImage, actualScan.getAnglesConstRef(), projectAlgo),
				  detWidth,
				  actualScan.getAnglesConstRef(),
				  0.0); //I0=0 because the Poisson statistics is not applied

		CTScan hi_dot = actualScan.getI0()*(-1*(actualScan/(actualScan.getI0()*(-1.0*li).exp())) + 1) * (-1.0*li).exp();
		//hi_dot.display("hi_dot");

		Phantom numerator("Numerator",
						  backProject(hi_dot, numberOfRecPoints, resolution, backProjectAlgo),
						  resolution);
		//numerator.display("hi_dot");

		//Calculate the curvatures
		CTScan curvatures("curvatures",
				(-1.0 * actualScan + actualScan.getI0()).getDataAsEigenMatrixRef().cwiseMax(0.0),
				detWidth,
				actualScan.getAnglesConstRef(),
				0.0); //I0=0 because the Poisson statistics is not applied
		//curvatures.display("Curvtures");

		Phantom denominator("Denominator",
							backProject(normFactors*curvatures, numberOfRecPoints, resolution, backProjectAlgo),
							resolution);
		//denominator.display();

		if(regularizerFunction == regularizerType::none){
			reconstImage = Phantom("reconstructedImage",
							        (reconstImage.getDataAsEigenMatrixRef().array() + (numerator.getDataAsEigenMatrixRef().array() )
			                              / (denominator.getDataAsEigenMatrixRef().array() ) ).cwiseMax(0.0),
										   resolution);
		}else if(regularizerFunction == regularizerType::quadratic){
			std::array<Phantom,2> regTerms = reconstImage.calculateQuadRegTerms();
			reconstImage = Phantom("reconstructedImage",
							               (reconstImage.getDataAsEigenMatrixRef().array() + (numerator.getDataAsEigenMatrixRef().array() - beta*regTerms[0].getDataAsEigenMatrixRef().array())
			                              / (denominator.getDataAsEigenMatrixRef().array() + beta*regTerms[1].getDataAsEigenMatrixRef().array()) ).cwiseMax(0.0),
										   resolution);
		}else if(regularizerFunction == regularizerType::Huber){
			for(int subIterNum=0; subIterNum < 1; ++subIterNum){
				//Calculate the regularization terms
				std::array<Phantom,2> regTerms = reconstImage.calculateHuberRegTerms(delta);

				reconstImage = Phantom("reconstructedImage",
				               (reconstImage.getDataAsEigenMatrixRef().array() + (numerator.getDataAsEigenMatrixRef().array() - beta*regTerms[0].getDataAsEigenMatrixRef().array())
                                                                               / (denominator.getDataAsEigenMatrixRef().array() + beta*regTerms[1].getDataAsEigenMatrixRef().array()) ).cwiseMax(0.0),
							   resolution);
			}
		}else if(regularizerFunction == regularizerType::Gibbs){
			for(int subIterNum=0; subIterNum < 1; ++subIterNum){
				//Calculate the regularization terms
				std::array<Phantom,2> regTerms = reconstImage.calculateGibbsRegTerms(delta);

				reconstImage = Phantom("reconstructedImage",
				                       (reconstImage.getDataAsEigenMatrixRef().array() + (numerator.getDataAsEigenMatrixRef().array() - beta*regTerms[0].getDataAsEigenMatrixRef().array())
			                                                                           / (denominator.getDataAsEigenMatrixRef().array() + beta*regTerms[1].getDataAsEigenMatrixRef().array()) ).cwiseMax(0.0),
				     				   resolution);
			}
		}
		else{
			std::cout << "\nERROR!! regularizer function not recognized, continue without regularization";
			reconstImage = Phantom("reconstructedImage",
										        (reconstImage.getDataAsEigenMatrixRef().array() + (numerator.getDataAsEigenMatrixRef().array() )
						                              / (denominator.getDataAsEigenMatrixRef().array() ) ).cwiseMax(0.0),
													   resolution);
		}

		if(referenceImage != ""){
			if(phantoms.find(referenceImage) == phantoms.end()){
				std::cout << std::endl << "ERROR!! referenceImage: \"" << referenceImage << "\" could not be found!! L2norm could not be calculated!";
				differenceNorms.push_back(0.0);
			} else{
				differenceNorms.push_back(reconstImage.compareNorm(phantoms[referenceImage]));
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "Reconstruction with SPS method took " << elapsedTime.count()/1000.0 << " seconds" << std::endl;

	//Move the backprojected image to reconsts map
	auto it = reconsts.find(imageID);
	if(it != reconsts.end()){
		std::cout << std::endl << "WARNING! A reconstructed image with label \"" << imageID << "\" already exists!!! Overwriting!!!";
		reconsts.erase(it);
	}
	reconsts.emplace(imageID, Reconst(imageID, reconstImage.getDataAsEigenMatrixRef(), resolution, differenceNorms));
}

#if ENABLE_CUDA
/***
 * Ray driven projection (Hao Gao method) implemented on the GPU
 * @param actualPhantom The phantom which should be projected
 * @param angles Vector of projection angles
 * @return Returns the sinogram of the actualPhantom
 */
Eigen::MatrixXd Gen1CT::project_rayDriven_GPU(const Phantom& actualPhantom,
		                    const Eigen::VectorXd& angles){
	std::cout << std::endl << "\nProjection with ray-driven method on GPU started";
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(static_cast<long>(pixNum), static_cast<long>(numAngles));

	auto pixSizes = actualPhantom.getPixSizes();
	auto numberOfPixels = actualPhantom.getNumberOfPixels();

	////////////////////////////////////////////////////////////////////////////////
	/// Projection with ray-driven method on GPU STARTS here !!
	////////////////////////////////////////////////////////////////////////////////

	launchRayDrivenProjectionKernel(actualPhantom.getDataAsEigenMatrixRef().data(), numberOfPixels, pixSizes,
			              numAngles, angles.data(), pixNum, detWidth,
						  sinogram.data());

	////////////////////////////////////////////////////////////////////////////////
	/// Projection with ray-driven method on GPU ENDS here !!
	////////////////////////////////////////////////////////////////////////////////

	auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "\nProjection with ray-driven method on GPU took " << elapsedTime.count() << " milliseconds\n";

	return sinogram;
}

Eigen::MatrixXd Gen1CT::backProject_rayDriven_GPU(const CTScan& sinogram,
						 const std::array<int,2>& numberOfRecPoints,
		                 const std::array<double,2>& resolution){

	std::cout << "\nBackprojection using ray-driven method on GPU started";
	auto start = std::chrono::high_resolution_clock::now();

	Eigen::MatrixXd backProjection = Eigen::MatrixXd::Zero(numberOfRecPoints[0], numberOfRecPoints[1]);

	const Eigen::VectorXd& angles=sinogram.getAnglesConstRef();
	int numAngles = angles.size();
	const Eigen::MatrixXd& sinoData=sinogram.getDataAsEigenMatrixRef();


	////////////////////////////////////////////////////////////////////////////////
	/// Backrojection with ray-driven method on GPU STARTS here !!
	////////////////////////////////////////////////////////////////////////////////

	launchRayDrivenBackprojectionKernel(sinoData.data(), numAngles, angles.data(), pixNum, detWidth, pixPositions.data(),
    										 backProjection.data(), numberOfRecPoints, resolution);

	////////////////////////////////////////////////////////////////////////////////
	/// Backrojection with ray-driven method on GPU ENDS here !!
	////////////////////////////////////////////////////////////////////////////////

	auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsedTime = stop - start;
	std::cout << "Backprojection using ray-driven method on GPU took " << elapsedTime.count() << " milliseconds\n";

	return backProjection;
}
#endif

double Gen1CT::compareReconToPhantom(std::string reconLabel, std::string phantomLabel) const{

	if(reconsts.find(reconLabel) == reconsts.end()){
						std::cout << std::endl << "ERROR!! reconLabel: \"" << reconLabel << "\" could not be found!! Abort mission";
						return 0.0;
	}

	if(phantoms.find(phantomLabel) == phantoms.end()){
						std::cout << std::endl << "ERROR!! phantomLabel: \"" << phantomLabel << "\" could not be found!! Abort mission";
						return 0.0;
	}

	return reconsts.at(reconLabel).compareNorm(phantoms.at(phantomLabel));
}

std::vector<double> Gen1CT::getConvergenceCurve(std::string label) const{
	std::vector<double> convergenceCurve(0);
	if(reconsts.find(label) == reconsts.end()){
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Can not return the convergence curve!";
	}else{
		convergenceCurve = reconsts.at(label).getConvergenceCurve();
	}

	return convergenceCurve;
}


