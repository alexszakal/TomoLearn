#include <TomoLearnS/Object2D.hpp>
#include <CImg.h>
#include <cstdint>
#include <iostream>
#include <vector>
#include <array>

Object2D::Object2D(const std::string& imageFilePath, const std::array<double, 2>& initObjPixSizes):
																objPixSizes{initObjPixSizes}{
	cimg_image=cimg_library::CImg<uint16_t>(imageFilePath.c_str());

	xPixCentresInMM = std::vector<double>( cimg_image._width  );
	yPixCentresInMM = std::vector<double>( cimg_image._height );

	objectSizeInMM = std::array<double, 2>{cimg_image._width * objPixSizes[0], cimg_image._height * objPixSizes[1]};

	std::cout << "Object loaded\n" << "size [mm*mm]: " << objectSizeInMM[0] << " * "<< objectSizeInMM[1]<< std::endl;

	for(unsigned int i=0; i<cimg_image._width; ++i){
		xPixCentresInMM[i]=-1*objectSizeInMM[0]/2 + i*objPixSizes[0] + objPixSizes[0]/2;
		yPixCentresInMM[i]=   objectSizeInMM[0]/2 - i*objPixSizes[1] - objPixSizes[1]/2;
	}
}

void Object2D::display(std::string title){
	cimg_window=cimg_library::CImgDisplay(cimg_image, title.c_str());
	//cimg_window.wait();
}

std::array<double,2> Object2D::getPixSizes(){
	return objPixSizes;
}

std::array<unsigned int, 2> Object2D::getNumberOfPixels(){
	return std::array<unsigned int, 2>{cimg_image._width, cimg_image._height};
}

double Object2D::linear_atY(int xPixelValue, double yCoordinateInMM){
	double yCoordinateInPixel = (objectSizeInMM[1]/2 - yCoordinateInMM) / objPixSizes[1];
	if (yCoordinateInPixel < 0.5){
		return 0;
	}
	if (yCoordinateInPixel > cimg_image._height-0.5) {
		return 0;
	}

	//Linear interpolation
	double neighbor0 = cimg_image(xPixelValue, floor(yCoordinateInPixel) );
	double neighbor1 = cimg_image(xPixelValue, ceil(yCoordinateInPixel) );

	return neighbor0 + (neighbor1 - neighbor0) * (yCoordinateInPixel - floor(yCoordinateInPixel));
}

double Object2D::linear_atX(int yPixelValue, double xCoordinateInMM){
	double xCoordinateInPixel = (objectSizeInMM[0]/2 + xCoordinateInMM) / objPixSizes[0];
	if (xCoordinateInPixel < 0.5){
			return 0;
	}
	if (xCoordinateInPixel > cimg_image._width-0.5) {
			return 0;
	}

	//Linear interpolation
	double neighbor0 = cimg_image(floor(xCoordinateInPixel), yPixelValue);
	double neighbor1 = cimg_image(ceil(xCoordinateInPixel), yPixelValue);

	return neighbor0 + (neighbor1 - neighbor0) * (xCoordinateInPixel - floor(xCoordinateInPixel));
}

double Object2D::getXValueAtPix(int pixValue){
	return xPixCentresInMM[pixValue];
}

double Object2D::getYValueAtPix(int pixValue){
	return yPixCentresInMM[pixValue];
}


