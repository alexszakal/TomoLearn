
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
			sinoImage(i,j) = static_cast<uint16_t>((sinogram(i,j)-minInt)/(maxInt-minInt)*65530);
		}
	}

	sinoWindow = cimg_library::CImgDisplay(sinoImage, "Sinogram");
	//sinoWindow.wait();

	//DEBUG
	/*matplotlibcpp::figure(44);
	Eigen::VectorXd proj = sinogram.col(0);
	matplotlibcpp::plot(std::vector<float> (&proj[0], proj.data()+proj.cols()*proj.rows()) );
	matplotlibcpp::show();
    */
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
	std::cout << std::endl << "Ratio of reconstructed and real (768. pixel): " << backprojection(768,821) / object->getDataAsEigenMatrixRef().row(821)(768);

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
	Eigen::VectorXd BPSlice = backprojection.col(821);
	Eigen::VectorXd ObjSlice = object->getDataAsEigenMatrixRef().row(821);
	matplotlibcpp::figure(27);
	matplotlibcpp::plot(std::vector<float> (&BPSlice[0], BPSlice.data()+BPSlice.cols()*BPSlice.rows()) );
	matplotlibcpp::plot(std::vector<float> (&ObjSlice[0], ObjSlice.data()+ObjSlice.cols()*ObjSlice.rows()) );
	matplotlibcpp::show(False);


}

void Gen1CT::FBP(std::vector<int> numberOfRecPoints, std::vector<double> resolution){
	/**
	 * Filtered BackProjection
	 * Input:
	 * 		-numberOfRecPoints [1]
	 * 		-resolution [mm]
	 * Output:
	 * 		None
	 */

	//Filtered Backprojection using Ram-Lak filter
	std::cout << "Filtering started" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	//Zero padding of the sinogram
	std::cout << "\n Eredeti pixnum: " << pixNum;
	int pixNumPadded = 2*pixNum-1;
	std::cout << "\n Padded pixnum: " << pixNumPadded;
	pixNumPadded = std::pow(2, std::ceil(std::log2(pixNumPadded)));
	std::cout << "\n Padded pixnum power2: " << pixNumPadded <<"\n";

	Eigen::MatrixXd paddedSinogram = Eigen::MatrixXd::Zero(pixNumPadded, numAngles);
	int startIndex=floor((pixNumPadded-pixNum)/2);
	paddedSinogram.block(startIndex, 0, pixNum, numAngles) = sinogram;

	//Fourier transform of the padded Sinogram:
	Eigen::FFT<double> fft;
	Eigen::MatrixXcd fftOfSinogram = Eigen::MatrixXcd::Zero(pixNumPadded, numAngles);
	for(int i=0; i<paddedSinogram.cols(); i++){
		fftOfSinogram.col(i) = fft.fwd(paddedSinogram.col(i));
	}

	//The Ram-Lak filter
	Eigen::ArrayXd freqFilter = Eigen::ArrayXd::Zero(pixNumPadded,1);
	for(int i=0; i<pixNumPadded; i++){
		if(pixNumPadded/2-std::abs(pixNumPadded/2-i) <= 700)
			freqFilter(i,0)=static_cast<double>( (pixNumPadded/2-std::abs(pixNumPadded/2-i)) )/(pixNumPadded/2);
		else
			freqFilter(i,0)=0;
	}
	//DEBUG plot the Ram-Lak
	matplotlibcpp::figure(3);
	matplotlibcpp::plot(std::vector<float> (&freqFilter[0], freqFilter.data()+freqFilter.cols()*freqFilter.rows()) );
	matplotlibcpp::show(False);


	//Multiply with filter
	for(int i=0; i<fftOfSinogram.cols(); ++i){
		fftOfSinogram.col(i) = fftOfSinogram.col(i).array() * freqFilter;
	}

	//IFFT of filtered sinogram
	for(int i=0; i<fftOfSinogram.cols(); i++){
		paddedSinogram.col(i) = fft.inv(fftOfSinogram.col(i)).real();
	}

	sinogram=detWidth/pixNum*paddedSinogram.block(startIndex,0, pixNum, numAngles);

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << "Filtering took " << duration.count() << " milliseconds" << std::endl;

}

const Eigen::MatrixXd& Gen1CT::getSinogram() const {
	return sinogram;
}

const std::vector<double>& Gen1CT::getPixPositions() const {
	return pixPositions;
}

void Gen1CT::FBP_bandlimited(std::vector<int> numberOfRecPoints, std::vector<double> resolution){
	/**
	 * Filtered BackProjection
	 * Input:
	 * 		-numberOfRecPoints [1]
	 * 		-resolution [mm]
	 * Output:
	 * 		None
	 */

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
	int maxFreq=256;
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

}

























