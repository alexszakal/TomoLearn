#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

#include <TomoLearnS/Gen1CT.hpp>
#include <TomoLearnS/Object2D.hpp>
#include <matplotlibcpp/matplotlibcpp.h>
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
	assert(object!=nullptr);

	angs=angles;
	numAngles = angles.size();

	sinogram = Eigen::MatrixXd::Zero(pixNum, numAngles);

	double t, theta;
	for(int pixI=0; pixI<pixNum; ++pixI){
		t=pixPositions[pixI];
		for(int angI=0; angI<numAngles; ++angI){
			theta = std::fmod(angles[angI], M_PI);
			double sinTheta = sin(theta);
			double cosTheta = cos(theta);

			if( (theta > M_PI/4 ) && (theta < 3*M_PI/4)  ){
				for(int objectXIndex=0; objectXIndex < object -> getNumberOfPixels()[0]; ++objectXIndex){
					double objectYinMM = t*sinTheta+ (t*cosTheta - object->getXValueAtPix(objectXIndex))/sinTheta*cosTheta;
					sinogram(pixI, angI) += object->linear_atY(objectXIndex, objectYinMM) / std::abs(sinTheta)*object->getPixSizes()[0];
				}
			} else{
				for(int objectYIndex=0; objectYIndex < object -> getNumberOfPixels()[1]; ++objectYIndex){
					double objectXinMM = t*cosTheta - (object->getYValueAtPix(objectYIndex)-t*sinTheta)/cosTheta*sinTheta;
					sinogram(pixI, angI) += object->linear_atX(objectYIndex, objectXinMM) / std::abs(cosTheta)*object->getPixSizes()[1];
				}
			}

		}
	}

}

void Gen1CT::displayMeasurement(){

	sinoImage=cimg_library::CImg<uint16_t>(pixNum, numAngles, 1, 1);

	double maxInt = sinogram.maxCoeff();
	for(int i=0; i<pixNum; ++i){
		for(int j=0; j<numAngles; ++j){
			sinoImage(i,j) = static_cast<uint16_t>(sinogram(i,j)/maxInt*65536);
		}
	}

	sinoWindow = cimg_library::CImgDisplay(sinoImage, "Sinogram");
	//sinoWindow.wait();
}

void Gen1CT::backProject(const std::vector<int>& numberOfRecPoints, const std::vector<double>& resolution){
	std::cout << "Reconstructrion started" << std::endl;
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
			for(unsigned int thIdx=0; thIdx<angs.size() ; ++thIdx){
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
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	std::cout << "Reconstruction took " << duration.count() << " seconds" << std::endl;

	BPImage = cimg_library::CImg<uint16_t>(numberOfRecPoints[0], numberOfRecPoints[1], 1, 1);
	double maxInt = backprojection.maxCoeff();
	for(int i=0; i<numberOfRecPoints[0]; ++i){
		for(int j=0; j<numberOfRecPoints[1]; ++j){
			BPImage(i,j) = static_cast<uint16_t>(backprojection(i,j)/maxInt*65536);
		}
	}

	BPWindow = cimg_library::CImgDisplay(BPImage, "Backprojection");

}

void Gen1CT::FBP(std::vector<int> numberOfRecPoints, std::vector<double> resolution){
	//Filtered Backprojection using Ram-Lak filter

	//Fourier transform of the Sinogram:
	Eigen::FFT<double> fft;
	Eigen::MatrixXcd fftOfSinogram = Eigen::MatrixXcd::Zero(pixNum, numAngles);
	for(int i=0; i<sinogram.cols(); i++){
		fftOfSinogram.col(i) = fft.fwd(sinogram.col(i));
		/*if(i==1){
			std::cout<<std::endl<<"Fourier transzformalt:\n";
			for(int j=0; j<pixNum; j++)
				std::cout<<fftOfSinogram(j,i)<<std::endl;
		}*/
	}

	//The Ram-Lak filter
	Eigen::ArrayXd freqFilter = Eigen::ArrayXd::Zero(pixNum,1);
	for(int i=0; i<pixNum; i++){
		if(pixNum/2-std::abs(pixNum/2-i) <= 1000)
			freqFilter(i,0)=pixNum/2-std::abs(pixNum/2-i);
		else
			freqFilter(i,0)=0;
	}
	//plot the Ram-Lak
	//matplotlibcpp::plot(std::vector<float> (&freqFilter[0], freqFilter.data()+freqFilter.cols()*freqFilter.rows()) );
	//matplotlibcpp::show();

	//Multiply with filter
	for(int i=0; i<fftOfSinogram.cols(); ++i){
		fftOfSinogram.col(i) = fftOfSinogram.col(i).array() * freqFilter;
	}

	//IFFT of filtered sinogram
	for(int i=0; i<fftOfSinogram.cols(); i++){
		sinogram.col(i) = fft.inv(fftOfSinogram.col(i));
	}



}



























