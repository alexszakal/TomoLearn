
#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

#include <TomoLearnS/Gen1CT.hpp>
#include <TomoLearnS/Object2D.hpp>
#include <matplotlibcpp/matplotlibcpp_old.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include <chrono>


Gen1CT::Gen1CT():detWidth{100},pixNum{100},object{nullptr},numAngles{0}{
	};

Gen1CT::Gen1CT(double detWidth, int pixNum):detWidth{detWidth},pixNum{pixNum},object{nullptr},numAngles{0}{
	pixPositions.resize(pixNum);
	double pixelDistance{detWidth/pixNum};
	for (int i=0; i<pixNum; i++){
		pixPositions[i]=-1*detWidth/2+(i+0.5)*pixelDistance;
	}
};

void Gen1CT::putObject(Object2D* sampleObject){
	object = sampleObject;
}

void Gen1CT::measure(const std::vector<double>& angles, int raysPerPixel){
	std::cout << "Radon transformation started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	assert(object!=nullptr);

	angs=angles;
	numAngles = angles.size();

	sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	std::array<double,2> pixSizes = object->getPixSizes();
	std::array<int,2> numberOfPixels = object -> getNumberOfPixels();

	double t, theta;
	const double piPer4 = M_PI/4;
	const double hpiPer4 = 3*M_PI/4;

	for(int pixI=0; pixI<pixNum; ++pixI){
		t=pixPositions[pixI];
		for(int angI=0; angI<numAngles; ++angI){
			theta = std::fmod(angles[angI], M_PI);
			double sinTheta = sin(theta);
			double cosTheta = cos(theta);
			double tanTheta = tan(theta);
			double cotTheta = 1/tanTheta;
			double absSinThetaInv = 1/std::abs(sinTheta);
			double absCosThetaInv = 1/std::abs(cosTheta);

			if( (theta > piPer4 ) && (theta < hpiPer4)  ){
				for(int objectXIndex=0; objectXIndex < numberOfPixels[0]; ++objectXIndex){
					double objectYinMM = t*sinTheta+ (t*cosTheta - object->getXValueAtPix(objectXIndex))*cotTheta;
					sinogram(pixI, angI) += object -> linear_atY(objectXIndex, objectYinMM) * absSinThetaInv*pixSizes[0];
				}
			} else{
				for(int objectYIndex=0; objectYIndex < numberOfPixels[1]; ++objectYIndex){
					double objectXinMM = t*cosTheta - (object->getYValueAtPix(objectYIndex)-t*sinTheta)*tanTheta;
					sinogram(pixI, angI) += object -> linear_atX(objectYIndex, objectXinMM) * absCosThetaInv*pixSizes[1];
				}
			}
		}
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
	std::cout << "Radon took " << duration.count() << " milliseconds" << std::endl;
}

void Gen1CT::displayMeasurement(){

	sinoImage=cimg_library::CImg<uint16_t>(pixNum, numAngles, 1, 1);

	double maxInt = sinogram.maxCoeff();
	double minInt = sinogram.minCoeff();
	for(int i=0; i<pixNum; ++i){
		for(int j=0; j<numAngles; ++j){
			sinoImage(i,j) = (static_cast<uint16_t>(sinogram(i,j)-minInt)/(maxInt-minInt)*65530);
		}
	}

	sinoWindow = cimg_library::CImgDisplay(sinoImage, "Sinogram");
	//sinoWindow.wait();
}

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
	for(int xIdx=0; xIdx<numberOfRecPoints[0]; ++xIdx){
		for(int yIdx=0; yIdx<numberOfRecPoints[1]; ++yIdx){
			//For every angle
			for(unsigned int thIdx=0; thIdx<angs.size(); ++thIdx){
				//Add the corresponding interpolated points from the sinogram
				double tValue = xValues[xIdx]*cosTheta[thIdx] + yValues[yIdx]*sinTheta[thIdx];
				if( (tValue<pixPositions[0]) || (tValue > pixPositions[pixNum-1]))
					continue;
				double pixIdx= (tValue + detWidth/2 - detWidth/pixNum/2) / (detWidth/pixNum);
				double valueInLowerPixel  = sinogram(floor(pixIdx), thIdx);
				double valueInHigherPixel = sinogram(ceil(pixIdx), thIdx);

				backprojection(xIdx, yIdx) += valueInLowerPixel + (valueInHigherPixel - valueInLowerPixel) * (pixIdx-floor(pixIdx));
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Reconstruction took " << duration.count() << " milliseconds" << std::endl;

	BPImage = cimg_library::CImg<uint16_t>(numberOfRecPoints[0], numberOfRecPoints[1], 1, 1);
	double maxInt = backprojection.maxCoeff();
	double minInt =backprojection.minCoeff();
	for(int i=0; i<numberOfRecPoints[0]; ++i){
		for(int j=0; j<numberOfRecPoints[1]; ++j){
			BPImage(i,j) = (static_cast<uint16_t>(backprojection(i,j)-minInt)/(maxInt-minInt)*65536);
		}
	}

	BPWindow = cimg_library::CImgDisplay(BPImage, "Backprojection");

	//DEBUG slice throughg the small ellipses
	Eigen::VectorXd slice = backprojection.col(821);
	matplotlibcpp::figure(2);
	matplotlibcpp::plot(std::vector<float> (&slice[0], slice.data()+slice.cols()*slice.rows()) );
	matplotlibcpp::show(False);

}

void Gen1CT::FBP(std::vector<int> numberOfRecPoints, std::vector<double> resolution){
	//Filtered Backprojection using Ram-Lak filter
	std::cout << "Filtering started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	//Fourier transform of the Sinogram:
	Eigen::FFT<double> fft;
	Eigen::MatrixXcd fftOfSinogram = Eigen::MatrixXcd::Zero(pixNum, numAngles);
	for(int i=0; i<sinogram.cols(); i++){
		fftOfSinogram.col(i) = fft.fwd(sinogram.col(i));
	}

	//The Ram-Lak filter
	//Eigen::ArrayXd freqFilter = Eigen::ArrayXd::Ones(pixNum,1)*3;
	Eigen::ArrayXd freqFilter = Eigen::ArrayXd::Zero(pixNum,1);
	for(int i=0; i<pixNum; i++){
		if(pixNum/2-std::abs(pixNum/2-i) <= 700)
			freqFilter(i,0)=static_cast<double>( (pixNum/2-std::abs(pixNum/2-i)) )/(pixNum/2) *2*M_PI/180;
		else
			freqFilter(i,0)=0;
	}
	/*//DEBUG plot the Ram-Lak
	matplotlibcpp::figure(3);
	matplotlibcpp::plot(std::vector<float> (&freqFilter[0], freqFilter.data()+freqFilter.cols()*freqFilter.rows()) );
	matplotlibcpp::show(False);
*/

	//Multiply with filter
	for(int i=0; i<fftOfSinogram.cols(); ++i){
		fftOfSinogram.col(i) = fftOfSinogram.col(i).array() * freqFilter;
	}

	//IFFT of filtered sinogram
	for(int i=0; i<fftOfSinogram.cols(); i++){
		sinogram.col(i) = fft.inv(fftOfSinogram.col(i)).real();
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Filtering took " << duration.count() << " milliseconds" << std::endl;

}



























