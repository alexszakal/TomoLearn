
#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

#include <Gen1CT.hpp>
#include <Phantom.hpp>
#include <Reconst.hpp>
#include <CTScan.hpp>
#include <Filter.hpp>

#include <matplotlibcpp.h>

#include <utility>   //std::swap
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <chrono>
#include <array>
#include <vector>
#include <random>
#include <algorithm> //std::min and std::max


#include <config.h>


    //REGI  pixPositionst le kell gyartani!!
Gen1CT::Gen1CT():detWidth{100},pixNum{100}{
	};


Gen1CT::Gen1CT(double detWidth, size_t pixNum):detWidth{detWidth},pixNum{pixNum}{
	pixPositions.resize( static_cast<size_t>(pixNum) );
	double pixelDistance{detWidth/pixNum};
	for (int i=0; i<pixNum; i++){
		pixPositions[ static_cast<size_t>(i) ]=-1*detWidth/2+(i+0.5)*pixelDistance;
	}
};

/**
 * Add a phantom from an existing image
 * @param label The label of the phantom in the library
 * @param phantomImageSource Location of the image file
 * @param pixSizes Pixel size in [mm]
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

void Gen1CT::displayPhantom(const std::string& label){
	if(phantoms.find(label) != phantoms.end()){
		phantoms[label].display(label);
	}
	else
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

void Gen1CT::setI0(double newI0){
	I0=newI0;
}

//void Gen1CT::measure_HaoGao(const std::string& phantomLabel,
//		                    const Eigen::VectorXd& angles,
//							const std::string& scanLabel){
//
//	if(phantoms.find(phantomLabel) == phantoms.end()){
//			std::cout << std::endl << "ERROR!! phantomLabel: \"" << phantomLabel << "\" could not be found!! Abort mission";
//			return;
//	}
//	Phantom& actualPhantom = phantoms[phantomLabel];
//
//	Eigen::MatrixXd sinogram = project_HaoGao_CPU(actualPhantom, angles);
//
//   	//Move the sinogram to CTScans map
//   	auto it = scans.find(scanLabel);
//   	if(it != scans.end()){
//   		std::cout << std::endl << "WARNING! A scan with label \"" << scanLabel << "\" already exists!!! Overwriting!!!";
//   		scans.erase(it);
//   	}
//   	scans.emplace(scanLabel, CTScan(scanLabel,sinogram, detWidth, angles));
//
//}


void Gen1CT::measure(const std::string& phantomLabel,
		     const Eigen::VectorXd& angles,
			 const std::string& scanLabel,
			 projectorType projector){

	if(phantoms.find(phantomLabel) == phantoms.end()){
		std::cout << std::endl << "ERROR!! phantomLabel: \"" << phantomLabel << "\" could not be found!! Abort mission";
		return;
	}

	Phantom& actualPhantom = phantoms[phantomLabel];
	int numAngles = angles.size();
	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	sinogram = project(actualPhantom, angles, projector);

	//Move the sinogram to CTScans map
	auto it = scans.find(scanLabel);
	if(it != scans.end()){
		std::cout << std::endl << "WARNING! A scan with label \"" << scanLabel << "\" already exists!!! Overwriting!!!";
		scans.erase(it);
	}

	Eigen::MatrixXd Isinogram;   //TODO Isinogram -> sinogram
	if(I0 !=0.0){ //Simulation of the effect of statistics on the line integrals has to be done
		//Calculate the expected value
		Isinogram = Eigen::exp( sinogram.array()* (-1.0) ) * I0 ;

		//Randomize the matrix
		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		std::mt19937_64 generator(seed1);    //Seed the generator
		for(int i=0; i<Isinogram.rows(); i++){
			for(int j=0; j<Isinogram.cols(); j++){
				std::poisson_distribution<int> poissonDist( Isinogram(i,j) );
				Isinogram(i,j) = (-1)* std::log( poissonDist(generator) / I0 );
				if (Isinogram(i,j)<0.0){    //The line integral has to be non-negative. It can become neg. because of statistics
					Isinogram(i,j) = 0.0;
				}
			}
		}
	}
	else{
		Isinogram = sinogram;
	}

	scans.emplace(scanLabel, CTScan(scanLabel,Isinogram, detWidth, angles, I0));
}

Eigen::MatrixXd Gen1CT::project(const Phantom& actualPhantom, const Eigen::VectorXd& angles, projectorType projector){
	int numAngles = angles.size();
	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

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
	}
	return sinogram;
}

Eigen::MatrixXd Gen1CT::project_rayDriven_CPU(const Phantom& actualPhantom,
		                    const Eigen::VectorXd& angles){
	std::cout << std::endl << "Projection with HaoGao's method started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(static_cast<long>(pixNum), static_cast<long>(numAngles));


	auto pixSizes = actualPhantom.getPixSizes();
    auto numberOfPixels = actualPhantom.getNumberOfPixels();

	//Convert Hounsfield to linear attenuation (LA) units
    Phantom actualPhantomLA = actualPhantom;                         //TODO: Torolni kell ha muxik randommal
//    if( I0 != 0.0){
//    	actualPhantomLA = (actualPhantom-1000.0) * (muWater/1000) + muWater;   // HU -> LA transform
//    	std::cout << "\n\nHounsfield units converted to linear attenuation units!!!";
//    }

    ////////////////////////////////////////////////////////////////////////////////
    /// HaoGao calculation STARTS here !!
    ////////////////////////////////////////////////////////////////////////////////

    //Coordinates of the column boundaries
    std::vector<double> colBoundaries( static_cast<unsigned long>(numberOfPixels[0]+1) );
    double halfPhantomWidth = numberOfPixels[0]*pixSizes[0]/2;
    double halfPhantomHeight = numberOfPixels[1]*pixSizes[1]/2;
    for(size_t boundIdx=0; boundIdx < static_cast<unsigned long>(numberOfPixels[0]+1); ++boundIdx){
    	colBoundaries[boundIdx] = -1*halfPhantomWidth+boundIdx*pixSizes[0];
    }

    //trigonometric functions for the angles
    std::vector<double>	thetaVector,
    	                sinThetaVector,
						cosThetaVector;

    for(int i=0; i<numAngles; i++){
    	thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
    	sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
    	cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
    }

    //Distance of the detector plane from origin which is outside of the phantom
    double detDist = 1.1 * sqrt(pow(pixSizes[0]*numberOfPixels[0], 2) + pow(pixSizes[1]*numberOfPixels[1], 2) ); //Distance between the detector plane and centre of rotation

    const double* dataPtr=actualPhantomLA.getDataAsEigenMatrixRef().data();

    double halfPixelSize=detWidth/pixNum/2;
    //Go through the angles
    for(size_t angI=0; angI < static_cast<size_t>(numAngles); ++angI){

    	//beam intersects the columns at most two pixels
    	if( pixSizes[1] / pixSizes[0] >= std::abs(std::tan(M_PI/2-thetaVector[angI])) ){

    		//go through the different t values
    		for(size_t pixIdx=0; pixIdx<pixNum; ++pixIdx){
    			double t = pixPositions[pixIdx] - halfPixelSize;
    			std::array<double,2> p1{ detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
    							                     -1*detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};
    			std::array<double,2> p2{ -1 * detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
    			                                     detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};

    			double ky = (p1[1]-p2[1])/(p1[0]-p2[0]);
    			//go through the columns of the image
    			for(int colIdx=0; colIdx<numberOfPixels[0]; ++colIdx){
    				double yi_minus = ( halfPhantomHeight - (ky*( colIdx   *pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) / pixSizes[1];
    				double yi_plus  = ( halfPhantomHeight - (ky*((colIdx+1)*pixSizes[0] - halfPhantomWidth - p1[0] ) + p1[1] ) ) / pixSizes[1];

    				//int Yi_minusIdx = static_cast<int>(yi_minus);
    				//int Yi_plusIdx = static_cast<int>(yi_plus);

    				int Yi_minusIdx = floor(yi_minus);
    				int Yi_plusIdx = floor(yi_plus);

    				double l, l_minus, l_plus; //intersecting length
    				if( Yi_minusIdx == Yi_plusIdx ){
    					if( (Yi_minusIdx < numberOfPixels[1]) and (Yi_minusIdx >= 0 ) ){
    						l=sqrt(1+ky*ky)*pixSizes[0];
//    						if( l * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx] < 0){  //DEBUG
//    							std::cout<<"\nHiba";
//    						}
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
    						//sinogram(pixIdx, angI) = l ; //DEBUG: why 0 when angI==45deg and pixI<256
    					}
    				} else{
    					if ( (Yi_minusIdx < numberOfPixels[1]) and (Yi_minusIdx >= 0 ) and (Yi_minusIdx < numberOfPixels[1]) ){
    						l_minus=(std::max(Yi_minusIdx, Yi_plusIdx)-yi_minus) / (yi_plus - yi_minus) * sqrt(1+ky*ky)*pixSizes[0];
//    						if( l_minus * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx] < 0){  //DEBUG
//    						    std::cout<<"\nHiba";
//    						}
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_minus * dataPtr[Yi_minusIdx*numberOfPixels[0] + colIdx];
    					}

    					if ( (Yi_plusIdx < numberOfPixels[1]) and (Yi_plusIdx >= 0 ) and (Yi_plusIdx < numberOfPixels[1]) ){
    						l_plus=(yi_plus - std::max(Yi_minusIdx, Yi_plusIdx)) / (yi_plus - yi_minus) * sqrt(1+ky*ky)*pixSizes[0];
//    						if( l_plus * dataPtr[Yi_plusIdx*numberOfPixels[0] + colIdx] < 0){  //DEBUG
//    							std::cout<<"\nHiba";
//    						}
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_plus * dataPtr[Yi_plusIdx*numberOfPixels[0] + colIdx];
    					}
    				}

    			}
    		}
    	}

    	//beam intersects the rows at most two pixels
    	if( pixSizes[1] / pixSizes[0] < std::abs(std::tan(M_PI/2-thetaVector[angI])) ){
    	    //go through the different t values
    	    for(size_t pixIdx=0; pixIdx<pixNum; ++pixIdx){
    	    	double t = pixPositions[pixIdx] - halfPixelSize;
    	    	std::array<double,2> p1{ detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
    	    					                     -1*detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};
    	    	std::array<double,2> p2{ -1 * detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
    	    	                                     detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};

    			double kx = (p1[0]-p2[0])/(p1[1]-p2[1]);
    			//go through the rows of the image
    	    	for(int rowIdx=0; rowIdx<numberOfPixels[1]; ++rowIdx){

    	    		double xi_minus = (halfPhantomWidth + (kx*( halfPhantomHeight - rowIdx     *pixSizes[1] - p1[1] ) + p1[0] ) ) / pixSizes[0];
    	    		double xi_plus  = (halfPhantomWidth + (kx*( halfPhantomHeight - (rowIdx+1) *pixSizes[1] - p1[1] ) + p1[0] ) ) / pixSizes[0];;

    	    		//int Xi_minusIdx = static_cast<int>(xi_minus); //statc_cast used instead of ceil to speed up calculation
    	    		//int Xi_plusIdx = static_cast<int>(xi_plus);   //statc_cast used instead of ceil to speed up calculation

    	    		int Xi_minusIdx = floor(xi_minus);
    	    		int Xi_plusIdx = floor(xi_plus);

    	    		double l, l_minus, l_plus; //intersecting length
    	    		if( Xi_minusIdx == Xi_plusIdx ){
    	    			if( (Xi_minusIdx < numberOfPixels[0]) and (Xi_minusIdx >= 0 ) ){
    	    				l=sqrt(1+kx*kx)*pixSizes[1];
//    	    				if( l * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx] < 0){  //DEBUG
//    	    					std::cout<<"\nHiba";
//    	    				}
    	    				sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
    	    				//sinogram(pixIdx, angI) = l ; //DEBUG: why 0 when angI==45deg and pixI<256
    	    			}
    	    		}
    	    		else{
    	    			if ( (Xi_minusIdx < numberOfPixels[0]) and (Xi_minusIdx >= 0 ) and (Xi_minusIdx < numberOfPixels[1]) ){
    	    				l_minus=(std::max(Xi_minusIdx, Xi_plusIdx)-xi_minus) / (xi_plus - xi_minus) * sqrt(1+kx*kx)*pixSizes[1];
//    	    				if( l_minus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx] < 0 ){  //DEBUG
//    	    				    std::cout<<"\nHiba";
//    	    				}
    	    				sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_minus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
    	    			}

    					if ( (Xi_plusIdx < numberOfPixels[0]) and (Xi_plusIdx >= 0 and (Xi_plusIdx < numberOfPixels[1])) ){
    						l_plus=(xi_plus - std::max(Xi_minusIdx, Xi_plusIdx)) / (xi_plus - xi_minus) * sqrt(1+kx*kx)*pixSizes[1];
//    						if( l_plus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx] < 0 ){  //DEBUG
//    						    std::cout<<"\nHiba";
//    						}
    						sinogram(static_cast<long>(pixIdx), static_cast<long>(angI)) += l_plus * dataPtr[rowIdx*numberOfPixels[0] + Xi_minusIdx];
    	    			}
    	    		}
    			}
    	    }
    	}
    }
    ////////////////////////////////////////////////////////////////////////////////
    /// HaoGao calculation ENDS here !!
    ////////////////////////////////////////////////////////////////////////////////

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
    std::cout << "Radon took " << duration.count() << " milliseconds" << std::endl;

//    Eigen::MatrixXd Isinogram;
//    //Simulate the counts with Poison statistics
//    if(I0 != 0.0){
//    	//Calculate the expected value
//    	Isinogram = Eigen::exp( sinogram.array()* (-1.0) ) * I0 ;
//
//    	//Randomize the matrix
//    	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
//    	std::mt19937_64 generator(seed1);    //Seed the generator
//    	for(int i=0; i<Isinogram.rows(); i++){
//    		for(int j=0; j<Isinogram.cols(); j++){
//    			std::poisson_distribution<int> poissonDist( Isinogram(i,j) );
//    			Isinogram(i,j) = (-1)* std::log( poissonDist(generator) / I0 );
//    		}
//    	}
//    }
//    else{
//    	Isinogram = sinogram;
//    }

    return sinogram;
}


Eigen::MatrixXd Gen1CT::project_pixelDriven_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles){

	std::cout << std::endl << "Projecting with interpolation started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	auto pixSizes = actualPhantom.getPixSizes();
	auto numberOfPixels = actualPhantom.getNumberOfPixels();

	//Convert Hounsfield to linear attenuation (LA) units
	Phantom actualPhantomLA = actualPhantom;                 //TODO: Torolni kell ha muxik randommal
//	if( I0 != 0.0){
//	  	actualPhantomLA = (actualPhantom-1000.0) * (muWater/1000) + muWater;   // HU -> LA transform
//	}

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
					double objectYinMM = t*sinThetaVector[angI]+ (t*cosThetaVector[angI] - actualPhantomLA.getXValueAtPix(objectXIndex))*cotThetaVector[angI];
					sinogram(static_cast<long>(pixI), static_cast<long>(angI) ) += actualPhantomLA.linear_atY(objectXIndex, objectYinMM) * absSinThetaInvVector[angI]*pixSizes[0];
				}
			} else{
				for(int objectYIndex=0; objectYIndex < numberOfPixels[1]; ++objectYIndex){
					double objectXinMM = t*cosThetaVector[angI] - (actualPhantomLA.getYValueAtPix(objectYIndex)-t*sinThetaVector[angI])*tanThetaVector[angI];
					sinogram(static_cast<long>(pixI), static_cast<long>(angI) ) += actualPhantomLA.linear_atX(objectYIndex, objectXinMM) * absCosThetaInvVector[angI]*pixSizes[1];
				}
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
	std::cout << "Projecting with interpolation took " << duration.count() << " milliseconds" << std::endl;

//	Eigen::MatrixXd Isinogram;
//	//Simulate the counts with Poison statistics
//	if(I0 != 0.0){
//		//Calculate the expected value
//		Isinogram = Eigen::exp( sinogram.array()* (-1.0) ) * I0 ;
//
//		//Randomize the matrix
//		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
//		std::mt19937_64 generator(seed1);    //Seed the generator
//		for(int i=0; i<Isinogram.rows(); i++){
//			for(int j=0; j<Isinogram.cols(); j++){
//				std::poisson_distribution<int> poissonDist( Isinogram(i,j) );
//				Isinogram(i,j) = (-1)* std::log( poissonDist(generator) / I0 );
//			}
//		}
//	}
//	else{
//		Isinogram = sinogram;
//	}

	return sinogram;

}

Eigen::MatrixXd Gen1CT::project_Siddon_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles){

	std::cout << std::endl << "Projection with improved Siddon's algo started";
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();
	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	auto pixSizes = actualPhantom.getPixSizes();
	auto numberOfPixels = actualPhantom.getNumberOfPixels();

	//Convert Hounsfield to linear attenuation (LA) units
	Phantom actualPhantomLA = actualPhantom;                      //TODO: Torolni kell ha muxik randommal
//	if( I0 != 0.0){
//		actualPhantomLA = (actualPhantom-1000.0) * (muWater/1000) + muWater;   // HU -> LA transform
//	}

	std::vector<double>	thetaVector,
						sinThetaVector, cosThetaVector;

	for(int i=0; i<numAngles; i++){
		thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
		sinThetaVector.push_back( sin(thetaVector[static_cast<size_t>(i)]) );
		cosThetaVector.push_back( cos(thetaVector[static_cast<size_t>(i)]) );
	}

	const double* dataPtr=actualPhantomLA.getDataAsEigenMatrixRef().data();

	double detDist = 1.1 * sqrt(pow(pixSizes[0]*numberOfPixels[0], 2) + pow(pixSizes[1]*numberOfPixels[1], 2) ); //Distance between the detector plane and centre of rotation

	double dX =actualPhantomLA.getPixSizes()[0];
	double dY =-1*actualPhantomLA.getPixSizes()[1];

	double bX = actualPhantomLA.getXValueAtPix(0) - dX/2;
	double bY = actualPhantomLA.getYValueAtPix(0) - dY/2;

	int Nx = actualPhantomLA.getNumberOfPixels()[0]+1;
	int Ny = actualPhantomLA.getNumberOfPixels()[1]+1;

	int iMin, iMax, jMin, jMax;

	//Improved Siddon's Algo:
	for(size_t angI=0; angI<numAngles; ++angI){
		for(size_t pixI=0; pixI<pixNum; ++pixI){
			double t=pixPositions[pixI];
			std::array<double,2> p1{ detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
									 -1*detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};
			std::array<double,2> p2{ -1 * detDist * sinThetaVector[angI] + t * cosThetaVector[angI],
									 detDist * cosThetaVector[angI] + t * sinThetaVector[angI]};

			double xDiff = p2[0] - p1[0];
			double yDiff = p2[1] - p1[1];

			double d_conv = sqrt(xDiff*xDiff + yDiff*yDiff);

			double d12 = 0; //Accumulator for the raySum

			if((angles[angI] == 0) || angles[angI] == M_PI){ //Edge case when ray is vertical
				int colIndex=std::floor( (p1[0]-bX)/dX );
				if (colIndex >=0 && colIndex < Nx-1){
					//Sum the matrix in vertical direction
					for(int rowIndex=0; rowIndex < Ny-1; ++rowIndex){
						d12 += dataPtr[ rowIndex*(Nx-1)+ colIndex];
					}
				}
				d12 *= -1*dY;
			}
			else if ((angles[angI] == M_PI/2) || angles[angI] == 3*M_PI/2){ //Edge case when ray is horizontal
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
						//double l_ij=(alpX-alp_c)*d_conv;
						//d12 += (alpX-alp_c) * dataPtr[ i + j*(Ny-1)];  Regi
						d12 += (alpX-alp_c) * dataPtr[ j*(Nx-1) + i];
						i += i_u;
						alp_c=alpX;
						alpX += alp_xu;
					} else{
						//double l_ij=(alpY-alp_c)*d_conv;
						//d12 += (alpY-alp_c) * dataPtr[ i + j*(Ny-1)];  Regi
						d12 += (alpY-alp_c) * dataPtr[ j*(Nx-1) + i];
						j += j_u;
						alp_c=alpY;
						alpY += alp_yu;
					}
				}
				d12 *= d_conv;
			}

			//Update the sinogram
			sinogram(pixI, angI) += d12;
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
	std::cout << "\nProjection with Siddon's algorithm took " << duration.count() << " milliseconds" << std::endl;


//	Eigen::MatrixXd Isinogram;
//	//Simulate the counts with Poison statistics
//	if(I0 != 0.0){
//		//Calculate the expected value
//		Isinogram = Eigen::exp( sinogram.array()* (-1.0) ) * I0 ;
//
//		//Randomize the matrix
//		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
//		std::mt19937_64 generator(seed1);    //Seed the generator
//		for(int i=0; i<Isinogram.rows(); i++){
//			for(int j=0; j<Isinogram.cols(); j++){
//				std::poisson_distribution<int> poissonDist( Isinogram(i,j) );
//				Isinogram(i,j) = (-1)* std::log( poissonDist(generator) / I0 );
//			}
//		}
//	}
//	else{
//		Isinogram = sinogram;
//	}

	return sinogram;
}

CTScan Gen1CT::getMeasurement(const std::string& label){
	if(scans.find(label) != scans.end()){
			return scans.at(label);
	}
	else{
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! returning an empty CTScan!!";
		return CTScan("EMPTY_SCAN", Eigen::MatrixXd::Zero(pixNum, 1), detWidth, Eigen::VectorXd{0}, I0);
	}
}

void Gen1CT::displayMeasurement(const std::string& label){
	if(scans.find(label) != scans.end()){
			scans.at(label).display(label);
		}
		else
			std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

Eigen::MatrixXd Gen1CT::backProject(const CTScan& sinogram,
									const std::array<int,2>& numberOfRecPoints,
									const std::array<double,2>& resolution,
									backprojectorType backProjector){

	switch(backProjector){
	case backprojectorType::pixelDriven:
		return backProject_pixelDriven_CPU(sinogram, numberOfRecPoints, resolution);
	case backprojectorType::rayDriven:
		return backProject_rayDriven_CPU(sinogram, numberOfRecPoints, resolution);
	}

}

/**
 * Backproject interpolated values of the sinogram to the pixels
 *
 * @param sinogram CTScan object which has to be backprojected
 * @param numberOfRecPoints Number of points in the backprojection
 * @param resolution Pixel size of the backprojected image
 * @return
 */
Eigen::MatrixXd Gen1CT::backProject_pixelDriven_CPU(const CTScan& sinogram,
						 const std::array<int,2>& numberOfRecPoints,
		                 const std::array<double,2>& resolution){

	std::cout << "Backprojection started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	assert(numberOfRecPoints.size()==2);
	assert(resolution.size()==2);

	Eigen::MatrixXd backprojection = Eigen::MatrixXd::Zero(numberOfRecPoints[0], numberOfRecPoints[1]);

	//Vectors with coordinates of the grid in real space
	double xMax=numberOfRecPoints[0]*resolution[0]/2;
	double xMin=-1*xMax;
	Eigen::VectorXd xValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[0], xMin, xMax)};

	double yMax=numberOfRecPoints[1]*resolution[1]/2;
	double yMin=-1*yMax;
	Eigen::VectorXd yValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[1], yMax, yMin)};

	const Eigen::VectorXd& angs=sinogram.getAnglesConstRef();
	const Eigen::MatrixXd& sinoData=sinogram.getDataAsEigenMatrixRef();

	std::vector<double> sinTheta(angs.size());
	std::vector<double> cosTheta(angs.size());
	for(unsigned int i=0; i<angs.size(); ++i){
		sinTheta[i]=sin(angs[i]);
		cosTheta[i]=cos(angs[i]);
	}

	//For each point in real space
		double offset = detWidth/2 + detWidth/pixNum/2; //BUG: itt +detWidth/pixNum/2 kell!!
		double invPixRes = pixNum/detWidth;
		double minPixPosition = -1*offset;
		double maxPixPosition = offset;
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
						if(pixIdxInt==0)
							valueInLowerPixel=0;
						else
							valueInLowerPixel  = sinoData(pixIdxInt-1, thIdx);
						valueInHigherPixel = sinoData(pixIdxInt, thIdx);
						backprojection(xIdx, yIdx) += valueInLowerPixel + (valueInHigherPixel - valueInLowerPixel) *
								                                           (1-(pixIdxInt-pixIdxDouble));
					}
					else{
						valueInLowerPixel  = sinoData(pixIdxInt, thIdx);
						if(pixIdxInt==pixNum-1)
							valueInHigherPixel = 0;
						else
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
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Reconstruction took " << duration.count() << " milliseconds" << std::endl;

	return backprojection;
}

/**
 * Backproject using the method proposed by Hao Gao in [REF]
 *
 * @param sinogram
 * @param numberOfRecPoints
 * @param resolution
 * @return
 */
Eigen::MatrixXd Gen1CT::backProject_rayDriven_CPU(const CTScan& sinogram,
						 const std::array<int,2>& numberOfRecPoints,
		                 const std::array<double,2>& resolution){

	std::cout << "Backprojection started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	assert(numberOfRecPoints.size()==2);
	assert(resolution.size()==2);

	Eigen::MatrixXd backProjection = Eigen::MatrixXd::Zero(numberOfRecPoints[0], numberOfRecPoints[1]);

	//Vectors with coordinates of the grid in real space
	double xMax=numberOfRecPoints[0]*resolution[0]/2;
	double xMin=-1*xMax;
	Eigen::VectorXd xValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[0], xMin, xMax)};

	double yMax=numberOfRecPoints[1]*resolution[1]/2;
	double yMin=-1*yMax;
	Eigen::VectorXd yValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[1], yMax, yMin)};

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
				//Determine the contributing detector pixels
				double xr = xValues(static_cast<long>(xIdx))*cosThetaVector[static_cast<size_t>(angIdx)]+
						       yValues(static_cast<long>(yIdx))*sinThetaVector[static_cast<size_t>(angIdx)];

				double minPixIdx = (xr + halfDetWidth - pixelRadius )*invDetPixSize;
				double maxPixIdx = (xr + halfDetWidth + pixelRadius )*invDetPixSize;

				//Calculate the intersection length and add to the backprojected image
				//Go through the possible pixels
				double lSum=0;
				double angleContrib = 0;
				for(int detPixIdx=std::max(0,static_cast<int>(minPixIdx));
					    detPixIdx < std::min(static_cast<int>(maxPixIdx),
							static_cast<int>(pixNum) ); ++detPixIdx){
					double a = cosThetaVector[static_cast<size_t>(angIdx)];
					double b = sinThetaVector[static_cast<size_t>(angIdx)];
					double c = pixPositions[static_cast<size_t>(detPixIdx)];

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
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Reconstruction took " << duration.count() << " milliseconds" << std::endl;

	return backProjection;
}

CTScan Gen1CT::applyFilter(const std::string& sinogramID, Filter filter){
	/**
	 * Filtering using Ram-Lak filter
	 * Input:
	 * 		-sinogramID
	 * Output:
	 * 		-CTScan object with the filtered sinogram
	 */

	//Timing
	std::cout << "Filtering started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	//Getting the reference to the sinogram to be filtered
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
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Filtering took " << duration.count() << " milliseconds" << std::endl;

	return CTScan(sinogramID+"F", paddedSinogram.block(startIndex,0, pixNum, numAngles),
		detWidth, sinogram.getAnglesConstRef(), sinogram.getI0()) ;
}

void Gen1CT::filteredBackProject(std::string sinogramID,
			                    const std::array<int,2>& numberOfRecPoints,
								const std::array<double,2>& resolution,
								FilterType filterType,
								double cutOffFreq,
								backprojectorType backProjectAlgo,
								const std::string& imageID){

	if(scans.find(sinogramID) == scans.end()){
			std::cout << std::endl << "ERROR!! sinogramID: \"" << sinogramID << "\" could not be found!! Abort mission";
			return;
	}

	//Apply filter on the sinogram
	CTScan filteredScan = applyFilter(sinogramID, Filter(filterType, cutOffFreq) );
	filteredScan.display(sinogramID + " filtered");

	//Backproject the image
	Eigen::MatrixXd backprojectedImage;
	backprojectedImage = backProject(filteredScan, numberOfRecPoints,resolution, backProjectAlgo);

	//Move the backprojected image to reconsts map
	auto it = reconsts.find(imageID);
	if(it != reconsts.end()){
		std::cout << std::endl << "WARNING! A reconstructed image with label \"" << imageID << "\" already exists!!! Overwriting!!!";
		reconsts.erase(it);
	}
	reconsts.emplace(imageID, Reconst(imageID, backprojectedImage, resolution));
}

void Gen1CT::displayReconstruction(const std::string& label){
	if(reconsts.find(label) != reconsts.end()){
		reconsts[label].display(label);
	}else
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

/**
 * Implement the MLEM algorithm
 */
void Gen1CT::MLEMReconst(std::string sinogramID,
			             const std::array<int,2>& numberOfRecPoints,
					     const std::array<double,2>& resolution,
						 projectorType projectAlgo,
						 backprojectorType backprojectAlgo,
						 const std::string& imageID,
						 int numberOfIterations){

	if(scans.find(sinogramID) == scans.end()){
				std::cout << std::endl << "ERROR!! sinogramID: \"" << sinogramID << "\" could not be found!! Abort mission";
				return;
	}
	CTScan actualScan = scans.at(sinogramID);

	//CT scan of const numbers:
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
	//normFactors.display("normFactors");

	//initialize the image
	Phantom reconstImage("reconstructedImage",
			             Eigen::MatrixXd::Ones(numberOfRecPoints[0], numberOfRecPoints[1])*muWater,
						 resolution);
//	reconstImage.display("ReconstImage");
//	std::cout << "Press ENTER to continue!";
//	std::cin.get();

	for(int itNumber = 0; itNumber < numberOfIterations; ++itNumber){

		std::cout << "\nIteration:" << itNumber+1 << " / " << numberOfIterations;

		//Forward project the estimation
		CTScan forwardProj( "iteration",
							project(reconstImage, actualScan.getAnglesConstRef(), projectAlgo),
							detWidth,
					        actualScan.getAnglesConstRef(),
							0.0); //I0=0 because the Poisson statistics is not applied

		if(forwardProj.getDataAsEigenMatrixRef().minCoeff() < 0){
			std::cout << "\n ERROR! negative value in forwardProj!";
		}

		//forwardProj.display("ForwardProj of Estimation");
		//std::cin.get();

		//Create the correction image
		//backproject the measurement elementwise divided with forwardprojection
		CTScan ratio = actualScan/(forwardProj+0.0005);
		CTScan tmpFwProj = forwardProj+0.0005;
//		tmpFwProj.display("tmpFwProj");
		//ratio.display("Ratio");
		//actualScan.display("actual Scan");

		Phantom corrImage;
		corrImage = Phantom("corrImage",
				           backProject( actualScan/(forwardProj+0.001),
								                   numberOfRecPoints,
				                                   resolution,
												   backprojectAlgo),
		                   resolution);
		//corrImage.display("Corrimage");   // ITT a corrImage(0,0) nagyon nagy.

		if(corrImage.getDataAsEigenMatrixRef().minCoeff() < 0){
					std::cout << "\n ERROR! negative value in corrImage!";
		}

		if(normFactors.getDataAsEigenMatrixRef().minCoeff() < 0){
							std::cout << "\n ERROR! negative value in normFactors!";
		}

		reconstImage = reconstImage / normFactors * corrImage;
		//reconstImage.display("reconstImage");

		//std::cout << "Iteration finished press ENTER to continue";
		//std::cin.get();

//		if(itNumber % 5 ==0)
//			reconstImage.display("Reconstruction state");
	}

	//reconstImage.display();

	//Move the backprojected image to reconsts map
	auto it = reconsts.find(imageID);
	if(it != reconsts.end()){
		std::cout << std::endl << "WARNING! A reconstructed image with label \"" << imageID << "\" already exists!!! Overwriting!!!";
		reconsts.erase(it);
	}
	reconsts.emplace(imageID, Reconst(imageID, reconstImage.getDataAsEigenMatrixRef(), resolution));
}

void Gen1CT::compareRowPhantomAndReconst(char direction, double position, const std::string& phantomID, const std::string& reconstID){
	/*** Compare the same rows of the phantom and the reconstruction
	 *
	 */
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

		std::vector<float> recXvals(recPixNum[1]);
		for(int idx=0; idx<recPixNum[1]; ++idx){
			recXvals[idx] = -1*recPixSizes[1]*recPixNum[1]/2 + (idx+0.5)*recPixSizes[1];
		}
		std::vector<float> phantomXvals(phantomPixNum[1]);
		for(int idx=0; idx<phantomPixNum[1]; ++idx){
			phantomXvals[idx] = -1*phantomPixSizes[1]*phantomPixNum[1]/2 + (idx+0.5)*phantomPixSizes[1];
		}

		matplotlibcpp::figure(27);
		matplotlibcpp::plot(recXvals, std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
		matplotlibcpp::plot(phantomXvals, std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
		matplotlibcpp::show();
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

		std::vector<float> recYvals(recPixNum[0]);
		for(int idx=0; idx<recPixNum[0]; ++idx){
			recYvals[idx] = recPixSizes[0]*recPixNum[0]/2 - (idx+0.5)*recPixSizes[0];
		}
		std::vector<float> phantomYvals(phantomPixNum[0]);
		for(int idx=0; idx<phantomPixNum[0]; ++idx){
			phantomYvals[idx] = phantomPixSizes[0]*phantomPixNum[0]/2 - (idx+0.5)*phantomPixSizes[0];
		}

		matplotlibcpp::figure(27);
		matplotlibcpp::plot(recYvals, std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
		matplotlibcpp::plot(phantomYvals, std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
		matplotlibcpp::show();
	}
	else{
		std::cout << "Unable to determine the slicing direction in compareRowPhantomAndReconst function."
				  << "\n direction should be \'X\' or \'Y\'. Aborting function";
		return;
	}
}

void Gen1CT::printPhantomParams(const std::string& phantomLabel){

	std::cout << "\n Phantom name: " << phantomLabel;
	std::cout << "\n Number of pixels: " << phantoms[phantomLabel].getNumberOfPixels()[0] << " x " << phantoms[phantomLabel].getNumberOfPixels()[1];
	std::cout << "\n Pixel sizes: " << phantoms[phantomLabel].getPixSizes()[0] << " x " << phantoms[phantomLabel].getPixSizes()[1];
}





















