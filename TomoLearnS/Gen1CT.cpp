
#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

#include <TomoLearnS/Gen1CT.hpp>
#include <TomoLearnS/Phantom.hpp>
#include <TomoLearnS/Reconst.hpp>
#include <TomoLearnS/CTScan.hpp>
#include <TomoLearnS/Filter.hpp>

#include <matplotlibcpp_old.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <chrono>
#include <array>
#include <vector>
#include <random>

    //REGI  pixPositionst le kell gyartani!!
Gen1CT::Gen1CT():detWidth{100},pixNum{100}{
	};


Gen1CT::Gen1CT(double detWidth, int pixNum):detWidth{detWidth},pixNum{pixNum}{
	pixPositions.resize( static_cast<size_t>(pixNum) );
	double pixelDistance{detWidth/pixNum};
	for (int i=0; i<pixNum; i++){
		pixPositions[ static_cast<size_t>(i) ]=-1*detWidth/2+(i+0.5)*pixelDistance;
	}
};

void Gen1CT::addPhantom(const std::string& label, const std::string& phantomImageSource, std::array<double, 2> pixSizes){
	/**
	 *   Add a phantom to the Phantom Library
	 */
	auto it = phantoms.find(label);
	if(it != phantoms.end()){
		std::cout << std::endl << "WARNING! A phantom with label \"" << label << "\"is already loaded!!! Overwriting!!!";
		phantoms.erase(it);
	}
	phantoms.emplace(label, Phantom(label,phantomImageSource,pixSizes));
}

void Gen1CT::addPhantom(const Phantom& newPhantom){
	/**
	 *   Add an existing Phantom to the class
	 */
	std::string phantomLabel = newPhantom.getLabel();
	auto it = phantoms.find(phantomLabel);
		if(it != phantoms.end()){
			std::cout << std::endl << "WARNING! A phantom with label \"" << phantomLabel << "\"is already loaded!!! Overwriting!!!";
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

void Gen1CT::measure_withInterpolation(const std::string& phantomLabel,
		             const Eigen::VectorXd& angles,
		             const std::string& scanLabel){

	std::cout << std::endl << "Radon transformation started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	if(phantoms.find(phantomLabel) == phantoms.end()){
		std::cout << std::endl << "ERROR!! phantomLabel: \"" << phantomLabel << "\" could not be found!! Abort mission";
		return;
	}

	Phantom& actualPhantom = phantoms[phantomLabel];

	auto pixSizes = actualPhantom.getPixSizes();
    auto numberOfPixels = actualPhantom.getNumberOfPixels();

	//Convert Hounsfield to linear attenuation (LA) units
    Phantom actualPhantomLA = actualPhantom;
    if( I0 != 0.0){
    	actualPhantomLA = (actualPhantom-1000.0) * (muWater/1000) + muWater;   // HU -> LA transform
    }

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
	std::cout << "Radon took " << duration.count() << " milliseconds" << std::endl;

	Eigen::MatrixXd Isinogram;
	//Simulate the counts with Poison statistics
	if(I0 != 0.0){
		//Calculate the expected value
		Isinogram = Eigen::exp( sinogram.array()* (-1.0) ) * I0 ;

		//Randomize the matrix
		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		std::mt19937_64 generator(seed1);    //Seed the generator
		for(int i=0; i<Isinogram.rows(); i++){
			for(int j=0; j<Isinogram.cols(); j++){
				std::poisson_distribution<int> poissonDist( Isinogram(i,j) );
				Isinogram(i,j) = (-1)* std::log( poissonDist(generator) / I0 );
			}
		}
	}
	else{
		Isinogram = sinogram;
	}

	//Move the sinogram to CTScans map
	auto it = scans.find(scanLabel);
	if(it != scans.end()){
		std::cout << std::endl << "WARNING! A scan with label \"" << phantomLabel << "\" already exists!!! Overwriting!!!";
		scans.erase(it);
	}
	scans.emplace(scanLabel, CTScan(scanLabel,Isinogram, detWidth, angles));
}

void Gen1CT::measure_Siddon(const std::string& phantomLabel,
		             const Eigen::VectorXd& angles,
		             const std::string& scanLabel){

	std::cout << std::endl << "Radon transformation with improved Siddon's algo started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

    int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	if(phantoms.find(phantomLabel) == phantoms.end()){
		std::cout << std::endl << "ERROR!! phantomLabel: \"" << phantomLabel << "\" could not be found!! Abort mission";
		return;
	}

	Phantom& actualPhantom = phantoms[phantomLabel];

	auto pixSizes = actualPhantom.getPixSizes();
    auto numberOfPixels = actualPhantom.getNumberOfPixels();

	//Convert Hounsfield to linear attenuation (LA) units
    Phantom actualPhantomLA = actualPhantom;
    if( I0 != 0.0){
    	actualPhantomLA = (actualPhantom-1000.0) * (muWater/1000) + muWater;   // HU -> LA transform
    }

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
	std::cout << "\n\nRadon with Siddon's algorithm took " << duration.count() << " milliseconds" << std::endl;


	Eigen::MatrixXd Isinogram;
	//Simulate the counts with Poison statistics
	if(I0 != 0.0){
		//Calculate the expected value
		Isinogram = Eigen::exp( sinogram.array()* (-1.0) ) * I0 ;

		//Randomize the matrix
		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		std::mt19937_64 generator(seed1);    //Seed the generator
		for(int i=0; i<Isinogram.rows(); i++){
			for(int j=0; j<Isinogram.cols(); j++){
				std::poisson_distribution<int> poissonDist( Isinogram(i,j) );
				Isinogram(i,j) = (-1)* std::log( poissonDist(generator) / I0 );
			}
		}
	}
	else{
		Isinogram = sinogram;
	}

	//Move the sinogram to CTScans map
	auto it = scans.find(scanLabel);
	if(it != scans.end()){
		std::cout << std::endl << "WARNING! A scan with label \"" << phantomLabel << "\" already exists!!! Overwriting!!!";
		scans.erase(it);
	}
	scans.emplace(scanLabel, CTScan(scanLabel,Isinogram, detWidth, angles));
}

CTScan Gen1CT::getMeasurement(const std::string& label){
	if(scans.find(label) != scans.end()){
			return scans.at(label);
	}
	else{
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! returning an empty CTScan!!";
		return CTScan("EMPTY_SCAN", Eigen::MatrixXd::Zero(pixNum, 1), detWidth, Eigen::VectorXd{0});
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
	double offset = detWidth/2 - detWidth/pixNum/2;
	double invPixRes = pixNum/detWidth;
	double minPixPosition = pixPositions[0];
	double maxPixPosition = pixPositions[pixNum-1];
	for(int xIdx=0; xIdx<numberOfRecPoints[0]; ++xIdx){
		for(int yIdx=0; yIdx<numberOfRecPoints[1]; ++yIdx){
			//For every angle
			for(unsigned int thIdx=0; thIdx<angs.size(); ++thIdx){
				//Add the corresponding interpolated points from the sinogram
				double tValue = xValues[xIdx]*cosTheta[thIdx] + yValues[yIdx]*sinTheta[thIdx];
				if( (tValue<minPixPosition) || (tValue > maxPixPosition))
					continue;
				double pixIdx= (tValue + offset) * invPixRes;
				double floorPixIdx = floor(pixIdx);
				double valueInLowerPixel  = sinoData(floorPixIdx, thIdx);
				double valueInHigherPixel = sinoData(ceil(pixIdx), thIdx);

				backprojection(xIdx, yIdx) += valueInLowerPixel + (valueInHigherPixel - valueInLowerPixel) * (pixIdx-floorPixIdx);
			}
		}
	}
	//Multiply with dTheta
	backprojection = backprojection*M_PI/angs.size();

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Reconstruction took " << duration.count() << " milliseconds" << std::endl;

	/* Regi kijelzes
	cimg_library::CImg<uint16_t> BPImage (numberOfRecPoints[0], numberOfRecPoints[1], 1, 1);
	double maxInt = backprojection.maxCoeff();
	double minInt =backprojection.minCoeff();
	for(int i=0; i<numberOfRecPoints[0]; ++i){
		for(int j=0; j<numberOfRecPoints[1]; ++j){
			BPImage(i,j) = static_cast<uint16_t>((backprojection(i,j)-minInt)/(maxInt-minInt)*65536);
		}
	}
	cimg_library::CImgDisplay BPImageDisplay;
	BPImage.display(BPImageDisplay, true,0, true);
	*/
	return backprojection;
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
		detWidth, sinogram.getAnglesConstRef()) ;
}

void Gen1CT::filteredBackProject(std::string sinogramID,
			                    const std::array<int,2>& numberOfRecPoints,
								const std::array<double,2>& resolution,
								FilterType filterType,
								double cutOffFreq,
								const std::string& imageID){

	if(scans.find(sinogramID) == scans.end()){
			std::cout << std::endl << "ERROR!! sinogramID: \"" << sinogramID << "\" could not be found!! Abort mission";
			return;
	}

	//Apply filter on the sinogram
	CTScan filteredScan = applyFilter(sinogramID, Filter(filterType, cutOffFreq) );
	filteredScan.display(sinogramID + " filtered");
	Eigen::MatrixXd backprojectedImage = backProject(filteredScan, numberOfRecPoints,resolution);

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

void Gen1CT::compareRowPhantomAndReconst(int rowNum, const std::string& phantomID, const std::string& reconstID){
	/*** Compare the same rows of the phantom and the reconstruction
	 *
	 */

	Eigen::VectorXd BPSlice = reconsts[reconstID].getDataAsEigenMatrixRef().col(rowNum);
	Eigen::VectorXd ObjSlice = phantoms[phantomID].getDataAsEigenMatrixRef().col(rowNum);
	matplotlibcpp::figure(27);
	matplotlibcpp::plot(std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
	matplotlibcpp::plot(std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
	matplotlibcpp::show();
}























