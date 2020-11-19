#include <TomoLearnS/Object2D.hpp>

#include <CImg.h>
#ifdef Success       //Because otherwise Eigen not compile
  #undef Success
#endif

#include <matplotlibcpp/matplotlibcpp_old.h>

#include <cstdint>
#include <iostream>
#include <vector>
#include <array>

Object2D::Object2D(const std::string& imageFilePath, const std::array<double, 2>& initObjPixSizes):
																objPixSizes{initObjPixSizes}{
	cimg_image = cimg_library::CImg<uint16_t>(imageFilePath.c_str());
	image = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(cimg_image._height, cimg_image._width);
	for(uint i=0; i<cimg_image._height; ++i){
		for(uint j=0; j<cimg_image._width; ++j){
			image(i,j) = cimg_image(j,i);   //Indexes are swapped due to different conventions in Eigen and cimg
		}
	}

	xPixCentresInMM = std::vector<double>( cimg_image._width  );
	yPixCentresInMM = std::vector<double>( cimg_image._height );

	objectSizeInMM = std::array<double, 2>{cimg_image._width * objPixSizes[0], cimg_image._height * objPixSizes[1]};

	std::cout << "Object loaded\n" << "size [mm*mm]: " << objectSizeInMM[0] << " * "<< objectSizeInMM[1]<< std::endl;

	for(unsigned int i=0; i<cimg_image._width; ++i){
		xPixCentresInMM[i]=-1*objectSizeInMM[0]/2 + i*objPixSizes[0] + objPixSizes[0]/2;
		yPixCentresInMM[i]=   objectSizeInMM[0]/2 - i*objPixSizes[1] - objPixSizes[1]/2;
	}

	/*//DEBUG show values along lower ellipses
	Eigen::VectorXd slice = image.row(821);
	matplotlibcpp::plot(std::vector<float> (&slice[0], slice.data()+slice.cols()*slice.rows()) );
	matplotlibcpp::show(False);
*/
}

void Object2D::display(std::string title){
	cimg_window=cimg_library::CImgDisplay(cimg_image, title.c_str());
	//cimg_window.wait();
}

std::array<double,2> Object2D::getPixSizes(){
	return objPixSizes;
}

std::array<int, 2> Object2D::getNumberOfPixels(){
	return std::array<int, 2>{static_cast<int>(cimg_image._width), static_cast<int>(cimg_image._height)};
}

double Object2D::linear_atY(int xPixelValue, double yCoordinateInMM){
	double yCoordinateInPixel = (objectSizeInMM[1]/2 - yCoordinateInMM) / objPixSizes[1];

	int lowerPixelIndex = floor(yCoordinateInPixel);
	int higherPixelIndex = ceil(yCoordinateInPixel);
	if( ( lowerPixelIndex < 0 ) || ( higherPixelIndex > static_cast<int>(cimg_image._height)-1 ) )
		return 0;

	//Linear interpolation
	double neighbor0 = cimg_image(xPixelValue, lowerPixelIndex );
	double neighbor1 = cimg_image(xPixelValue, higherPixelIndex );

	return neighbor0 + (neighbor1 - neighbor0) * (yCoordinateInPixel - lowerPixelIndex);
}

double Object2D::linear_atX(int yPixelValue, double xCoordinateInMM){
	double xCoordinateInPixel = (objectSizeInMM[0]/2 + xCoordinateInMM) / objPixSizes[0];
	int lowerPixelIndex = floor(xCoordinateInPixel);
	int higherPixelIndex = ceil(xCoordinateInPixel);
	if( ( lowerPixelIndex < 0 ) || ( higherPixelIndex > static_cast<int>(cimg_image._height)-1 ) )
		return 0;

	//Linear interpolation
	double neighbor0 = cimg_image(lowerPixelIndex, yPixelValue);
	double neighbor1 = cimg_image(higherPixelIndex, yPixelValue);

	return neighbor0 + (neighbor1 - neighbor0) * (xCoordinateInPixel - lowerPixelIndex);
}

double Object2D::getXValueAtPix(int pixValue){
	return xPixCentresInMM[pixValue];
}

double Object2D::getYValueAtPix(int pixValue){
	return yPixCentresInMM[pixValue];
}

const Eigen::MatrixXd& Object2D::getDataAsEigenMatrixRef() const{
	return image;
}


