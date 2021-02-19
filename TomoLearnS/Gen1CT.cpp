
#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

#include <TomoLearnS/Gen1CT.hpp>
//#include <TomoLearnS/Object2D.hpp>
#include <TomoLearnS/Phantom.hpp>
#include <TomoLearnS/Reconst.hpp>
#include <TomoLearnS/CTScan.hpp>

#include <matplotlibcpp/matplotlibcpp_old.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include <chrono>
#include <array>
#include <vector>

    //REGI  pixPositionst le kell gyartani!!
Gen1CT::Gen1CT():detWidth{100},pixNum{100}{
	};


Gen1CT::Gen1CT(double detWidth, int pixNum):detWidth{detWidth},pixNum{pixNum}{
	pixPositions.resize(pixNum);
	double pixelDistance{detWidth/pixNum};
	for (int i=0; i<pixNum; i++){
		pixPositions[i]=-1*detWidth/2+(i+0.5)*pixelDistance;
	}
};

void Gen1CT::addPhantom(const std::string& label, const std::string& phantomImageSource, std::array<double, 2> pixSizes){
	/** Add a phantom to the Phantom Library
	 *
	 */
	auto it = phantoms.find(label);
	if(it != phantoms.end()){
		std::cout << std::endl << "WARNING! A phantom with label \"" << label << "\"is already loaded!!! Overwriting!!!";
		phantoms.erase(it);
	}
	phantoms.emplace(label, Phantom(label,phantomImageSource,pixSizes));
}


void Gen1CT::displayPhantom(const std::string& label, const std::string& title){
	if(phantoms.find(label) != phantoms.end()){
		phantoms.at(label).display();
	}
	else
		std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

void Gen1CT::measure(const std::string& phantomLabel,
		             const Eigen::VectorXd& angles,
		             const std::string& scanLabel){
	/* Idomeres kezdete helyesen
	std::cout << std::endl << "Radon transformation started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
*/

	int numAngles = angles.size();

	Eigen::MatrixXd sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	if(phantoms.find(phantomLabel) == phantoms.end()){
		std::cout << std::endl << "ERROR!! phantomLabel: \"" << phantomLabel << "\" could not be found!! Abort mission";
		return;
	}
	Phantom& actualPhantom = phantoms.at(phantomLabel);

	auto pixSizes = phantoms.at(phantomLabel).getPixSizes();
	auto numberOfPixels = phantoms.at(phantomLabel).getNumberOfPixels();

	double t;
	const double piPer4 = M_PI/4;
	const double hpiPer4 = 3*M_PI/4;

	std::vector<double>	thetaVector, sinThetaVector, cosThetaVector, tanThetaVector, cotThetaVector,
	                    absSinThetaInvVector, absCosThetaInvVector;
	for(int i=0; i<numAngles; i++){
		thetaVector.push_back( std::fmod(angles[i], 2*M_PI) );
		sinThetaVector.push_back( sin(thetaVector[i]) );
		cosThetaVector.push_back( cos(thetaVector[i]) );
		tanThetaVector.push_back( tan(thetaVector[i]) );
		cotThetaVector.push_back( 1/tanThetaVector[i] );
		absSinThetaInvVector.push_back( 1/std::abs(sinThetaVector[i]) );
		absCosThetaInvVector.push_back( 1/std::abs(cosThetaVector[i]) );
	}

//Idomeres ideiglenes kezdete
std::cout << std::endl << "Radon transformation started" << std::endl;
auto start = std::chrono::high_resolution_clock::now();
	for(int pixI=0; pixI<pixNum; ++pixI){

		t=pixPositions[pixI];
		for(int angI=0; angI<numAngles; ++angI){

			if( (thetaVector[angI] > piPer4 ) && (thetaVector[angI] < hpiPer4)  ){
				for(int objectXIndex=0; objectXIndex < numberOfPixels[0]; ++objectXIndex){
					double objectYinMM = t*sinThetaVector[angI]+ (t*cosThetaVector[angI] - actualPhantom.getXValueAtPix(objectXIndex))*cotThetaVector[angI];
					sinogram(pixI, angI) += actualPhantom.linear_atY(objectXIndex, objectYinMM) * absSinThetaInvVector[angI]*pixSizes[0];
				}
			} else{
				for(int objectYIndex=0; objectYIndex < numberOfPixels[1]; ++objectYIndex){
					double objectXinMM = t*cosThetaVector[angI] - (actualPhantom.getYValueAtPix(objectYIndex)-t*sinThetaVector[angI])*tanThetaVector[angI];
					sinogram(pixI, angI) += actualPhantom.linear_atX(objectYIndex, objectXinMM) * absCosThetaInvVector[angI]*pixSizes[1];
				}
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
	std::cout << "Radon took " << duration.count() << " milliseconds" << std::endl;

	//Move the sinogram to CTScans map
	auto it = scans.find(scanLabel);
	if(it != scans.end()){
		std::cout << std::endl << "WARNING! A scan with label \"" << phantomLabel << "\" already exists!!! Overwriting!!!";
		scans.erase(it);
	}
	scans.emplace(scanLabel, CTScan(scanLabel,sinogram, detWidth, angles));
}


/*   Mukodo Radon iteracio
 * for(int pixI=0; pixI<pixNum; ++pixI){

		t=pixPositions[pixI];
		for(int angI=0; angI<numAngles; ++angI){

			if( (thetaVector[angI] > piPer4 ) && (thetaVector[angI] < hpiPer4)  ){
				for(int objectXIndex=0; objectXIndex < numberOfPixels[0]; ++objectXIndex){
					double objectYinMM = t*sinThetaVector[angI]+ (t*cosThetaVector[angI] - phantoms.at(phantomLabel).getXValueAtPix(objectXIndex))*cotThetaVector[angI];
					sinogram(pixI, angI) += phantoms.at(phantomLabel).linear_atY(objectXIndex, objectYinMM) * absSinThetaInvVector[angI]*pixSizes[0];
				}
			} else{
				for(int objectYIndex=0; objectYIndex < numberOfPixels[1]; ++objectYIndex){
					double objectXinMM = t*cosThetaVector[angI] - (phantoms.at(phantomLabel).getYValueAtPix(objectYIndex)-t*sinThetaVector[angI])*tanThetaVector[angI];
					sinogram(pixI, angI) += phantoms.at(phantomLabel).linear_atX(objectYIndex, objectXinMM) * absCosThetaInvVector[angI]*pixSizes[1];
				}
			}
		}
	}
 */







void Gen1CT::displayMeasurement(const std::string& label){
	if(scans.find(label) != scans.end()){
			scans.at(label).display();
		}
		else
			std::cout << std::endl << "ERROR!! Label: \"" << label << "\" could not be found!! Skipping the display.";
}

//REGI
/*
void Gen1CT::backProject(const std::vector<int>& numberOfRecPoints, const std::vector<double>& resolution){
	std::cout << "Backprojection started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	assert(numberOfRecPoints.size()==2);
	assert(resolution.size()==2);

	backprojection = Eigen::MatrixXd::Zero(numberOfRecPoints[0], numberOfRecPoints[1]);

	//Vectors with coordinates of the grid in real space
	double xMax=numberOfRecPoints[0]*resolution[0]/2;
	double xMin=-1*xMax;
	Eigen::VectorXd xValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[0], xMin, xMax)};

	double yMax=numberOfRecPoints[1]*resolution[1]/2;
	double yMin=-1*yMax;
	Eigen::VectorXd yValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[1], yMax, yMin)};

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
				double valueInLowerPixel  = sinogram(floorPixIdx, thIdx);
				double valueInHigherPixel = sinogram(ceil(pixIdx), thIdx);

				backprojection(xIdx, yIdx) += valueInLowerPixel + (valueInHigherPixel - valueInLowerPixel) * (pixIdx-floorPixIdx);
			}
		}
	}
	//Multiply with dTheta
	backprojection = backprojection*M_PI/angs.size();

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Reconstruction took " << duration.count() << " milliseconds" << std::endl;

	//DEBUG
	//for( int i=0; i<numberOfRecPoints[1]; ++i){
	//		std::cout << std::endl << "Index: "<< i << "  Image: " << object->getDataAsEigenMatrixRef().row(821)(i)   << "  Reconst: " << backprojection(i,821);
	//}

	//TODO mukodjon a kovetkezo sor
	//std::cout << std::endl << "Ratio of reconstructed and real (768. pixel): " << backprojection(768,821) / object->getDataAsEigenMatrixRef().row(821)(768);

	BPImage = cimg_library::CImg<uint16_t>(numberOfRecPoints[0], numberOfRecPoints[1], 1, 1);
	double maxInt = backprojection.maxCoeff();
	double minInt =backprojection.minCoeff();
	for(int i=0; i<numberOfRecPoints[0]; ++i){
		for(int j=0; j<numberOfRecPoints[1]; ++j){
			BPImage(i,j) = static_cast<uint16_t>((backprojection(i,j)-minInt)/(maxInt-minInt)*65536);
		}
	}

	BPWindow = cimg_library::CImgDisplay(BPImage, "Backprojection");

	//DEBUG slice through the small ellipses
	//TODO Mukodjon ez a teszt. Inkabb kulon legyen a main()-ben
  //Eigen::VectorXd BPSlice = backprojection.col(821);
  //Eigen::VectorXd ObjSlice = object->getDataAsEigenMatrixRef().row(821);
  //matplotlibcpp::figure(27);
  //matplotlibcpp::plot(std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
  //matplotlibcpp::plot(std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
  //matplotlibcpp::show(False);
}
*/


void Gen1CT::FBP(std::vector<int> numberOfRecPoints, std::vector<double> resolution){
	/**
	 * Filtered BackProjection
	 * Input:
	 * 		-numberOfRecPoints [1]
	 * 		-resolution [mm]
	 * Output:
	 * 		None
	 */
/*
	//Filtered Backprojection using Ram-Lak filter
	//Timing
	std::cout << "Filtering started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	//Zero padding of the sinogram
	int pixNumPadded = std::pow(2, std::ceil(std::log2(2*pixNum-1)));
	std::cout << "\nOriginal pixNum:" << pixNum << " ZeroPadded to: " << pixNumPadded <<"\n";

	Eigen::MatrixXd paddedSinogram = Eigen::MatrixXd::Zero(pixNumPadded, numAngles);
	int startIndex=floor((pixNumPadded-pixNum)/2);
	paddedSinogram.block(startIndex, 0, pixNum, numAngles) = sinogram;

	//DEBUG
	cimg_library::CImg<uint16_t> paddedImage = cimg_library::CImg<uint16_t>(pixNumPadded, numAngles, 1, 1);
	double maxInt = paddedSinogram.maxCoeff();
	double minInt = paddedSinogram.minCoeff();
	for(int i=0; i<pixNumPadded; ++i){
		for(int j=0; j<numAngles; ++j){
			paddedImage(i,j) = static_cast<uint16_t>((paddedSinogram(i,j)-minInt)/(maxInt-minInt)*65530);
		}
	}
	cimg_library::CImgDisplay paddedImageDisplay = cimg_library::CImgDisplay(paddedImage, "PaddedImage");
	paddedImageDisplay.wait();
	//DEBUG END

	//Fourier transform of the padded Sinogram:
	Eigen::FFT<double> fft;
	Eigen::MatrixXcd fftOfSinogram = Eigen::MatrixXcd::Zero(pixNumPadded, numAngles);
	for(int i=0; i<paddedSinogram.cols(); i++){
		fftOfSinogram.col(i) = fft.fwd(paddedSinogram.col(i));
	}

	//Construct the filter
	Eigen::MatrixXcd freqFilter = Eigen::MatrixXd::Zero(pixNumPadded,1);
	Eigen::MatrixXd filter = Eigen::MatrixXd::Zero(pixNumPadded,1);

	double tau=detWidth/pixNum;
	filter(0)=1/(4*tau*tau);
	for(int i=0; i<pixNumPadded/4; ++i){
		int idx= i*2+1;
		filter(idx) = -1/std::pow(idx*M_PI*tau, 2);
		filter(pixNumPadded-idx) = filter(idx);
	}
	Eigen::FFT<double> fft2;
	freqFilter.col(0)=fft2.fwd(filter.col(0));

	//Low-pass filter to suppress the noise
	int maxFreq=1024;
	for(int i=maxFreq+1; i<=pixNumPadded/2; ++i ){
		freqFilter(i)=0;
		freqFilter(pixNumPadded-i)=0;
	}

	//DEBUG plot the Ram-Lak
	matplotlibcpp::figure(3);
	Eigen::MatrixXd absVector = freqFilter.imag().array().pow(2) + freqFilter.real().array().pow(2) ;
	absVector = absVector.array().pow(0.5);
	std::cout << "\n H(0)= " << absVector(0);
	matplotlibcpp::plot(std::vector<float> (&absVector(0), absVector.data()+absVector.cols()*absVector.rows()) );
	matplotlibcpp::show();


	//Multiply with filter
	for(int i=0; i<fftOfSinogram.cols(); ++i){
		fftOfSinogram.col(i) = fftOfSinogram.col(i).array() * freqFilter.array();
	}

	//IFFT of filtered sinogram
	for(int i=0; i<fftOfSinogram.cols(); i++){
		paddedSinogram.col(i) = tau * fft.inv(fftOfSinogram.col(i)).real();
	}

	sinogram=paddedSinogram.block(startIndex,0, pixNum, numAngles);

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Filtering took " << duration.count() << " milliseconds" << std::endl;

	std::cout << "Backprojection started" << std::endl;
		auto start = std::chrono::high_resolution_clock::now();

		assert(numberOfRecPoints.size()==2);
		assert(resolution.size()==2);

		backprojection = Eigen::MatrixXd::Zero(numberOfRecPoints[0], numberOfRecPoints[1]);

		//Vectors with coordinates of the grid in real space
		double xMax=numberOfRecPoints[0]*resolution[0]/2;
		double xMin=-1*xMax;
		Eigen::VectorXd xValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[0], xMin, xMax)};

		double yMax=numberOfRecPoints[1]*resolution[1]/2;
		double yMin=-1*yMax;
		Eigen::VectorXd yValues{Eigen::VectorXd::LinSpaced(numberOfRecPoints[1], yMax, yMin)};

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
					double valueInLowerPixel  = sinogram(floorPixIdx, thIdx);
					double valueInHigherPixel = sinogram(ceil(pixIdx), thIdx);

					backprojection(xIdx, yIdx) += valueInLowerPixel + (valueInHigherPixel - valueInLowerPixel) * (pixIdx-floorPixIdx);
				}
			}
		}
		//Multiply with dTheta
		backprojection = backprojection*M_PI/angs.size();

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
		std::cout << "Reconstruction took " << duration.count() << " milliseconds" << std::endl;

		//DEBUG
		//for( int i=0; i<numberOfRecPoints[1]; ++i){
		//		std::cout << std::endl << "Index: "<< i << "  Image: " << object->getDataAsEigenMatrixRef().row(821)(i)   << "  Reconst: " << backprojection(i,821);
		//}

		//TODO mukodjon a kovetkezo sor
		//std::cout << std::endl << "Ratio of reconstructed and real (768. pixel): " << backprojection(768,821) / object->getDataAsEigenMatrixRef().row(821)(768);

		BPImage = cimg_library::CImg<uint16_t>(numberOfRecPoints[0], numberOfRecPoints[1], 1, 1);
		double maxInt = backprojection.maxCoeff();
		double minInt =backprojection.minCoeff();
		for(int i=0; i<numberOfRecPoints[0]; ++i){
			for(int j=0; j<numberOfRecPoints[1]; ++j){
				BPImage(i,j) = static_cast<uint16_t>((backprojection(i,j)-minInt)/(maxInt-minInt)*65536);
			}
		}

		BPWindow = cimg_library::CImgDisplay(BPImage, "Backprojection");

		//DEBUG slice through the small ellipses
		//TODO Mukodjon ez a teszt. Inkabb kulon legyen a main()-ben
//		Eigen::VectorXd BPSlice = backprojection.col(821);
//		Eigen::VectorXd ObjSlice = object->getDataAsEigenMatrixRef().row(821);
//		matplotlibcpp::figure(27);
//		matplotlibcpp::plot(std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
//		matplotlibcpp::plot(std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
//		matplotlibcpp::show(False);

*/
}

























