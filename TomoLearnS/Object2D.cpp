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

Object2D::Object2D(const std::string& label, const std::string& imageFilePath, const std::array<double, 2>& objPixSizes):
																	objPixSizes{objPixSizes}, label{label}{
	//Read the image file
    cimg_library::CImg<uint16_t> imageLoader;
	imageLoader = cimg_library::CImg<uint16_t>(imageFilePath.c_str());
	//Copy data to Eigen::MatrixXd
	objData = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(imageLoader._width, imageLoader._height);
	for(uint i=0; i<imageLoader._width; ++i){
		for(uint j=0; j<imageLoader._height; ++j){
			objData(i,j) = imageLoader(i,j);
		}
	}

	numberOfPixels = {static_cast<int>(objData.rows()), static_cast<int>(objData.cols())};
	xPixCentreCoords = std::vector<double>( objData.cols()  );
	yPixCentreCoords = std::vector<double>( objData.rows() );

	objWidthHeight = std::array<double, 2>{numberOfPixels[0] * objPixSizes[0], numberOfPixels[1] * objPixSizes[1]};

	std::cout << std:: endl << "Object in \"" << imageFilePath << "\" loaded" << std::endl <<
			"size [mm*mm]: " << objWidthHeight[0] << " * "<< objWidthHeight[1]<< std::endl;

	for(unsigned int i=0; i<imageLoader._width; ++i){
		xPixCentreCoords[i]=-1*objWidthHeight[0]/2 + i*objPixSizes[0] + objPixSizes[0]/2;
		yPixCentreCoords[i]=   objWidthHeight[0]/2 - i*objPixSizes[1] - objPixSizes[1]/2;
	}

	/*//DEBUG show values along lower ellipses
	Eigen::VectorXd slice = image.row(821);
	matplotlibcpp::plot(std::vector<float> (&slice[0], slice.data()+slice.cols()*slice.rows()) );
	matplotlibcpp::show(False);
*/
}

void Object2D::display(std::string title){
	std::stringstream ss;
	ss << title << " Label: " << label << " " << numberOfPixels[0] << " x " << numberOfPixels[1] << " points; "
			                                     << objPixSizes[0] << " x " << objPixSizes[1] << " pixSize"
									             << objWidthHeight[0] << " x " << objWidthHeight[1] << " width x height";
	title = ss.str();

	cimg_image=cimg_library::CImg<uint16_t>(numberOfPixels[0], numberOfPixels[1], 1, 1);

	double maxInt = objData.maxCoeff();
	double minInt = objData.minCoeff();
	double normFactor = 65530/(maxInt-minInt);
	for(int i=0; i<numberOfPixels[0]; ++i){
		for(int j=0; j<numberOfPixels[1]; ++j){
			cimg_image(i,j) = static_cast<uint16_t>((objData(i,j)-minInt)*normFactor);
		}
	}

	cimg_window = cimg_library::CImgDisplay(cimg_image, title.c_str());
}

std::array<double,2> Object2D::getPixSizes() const{
	/** Return the size of pixels in units of the coordinate system
	 *
	 */
	return objPixSizes;
}

std::array<int, 2> Object2D::getNumberOfPixels() const{
	/** Returns the number of pixels
	 *
	 */
	return numberOfPixels;
}

double Object2D::linear_atY(int xPixelValue, double yCoordinateInMM) const{
	/** Linear interpolation in X direction at yCoordinateInMM exactly at xPixelValue
	 *
	 */
	double yCoordinateInPixel = (objWidthHeight[1]/2 - yCoordinateInMM) / objPixSizes[1];

	int lowerPixelIndex = floor(yCoordinateInPixel);
	int higherPixelIndex = ceil(yCoordinateInPixel);
	if( ( lowerPixelIndex < 0 ) || ( higherPixelIndex > static_cast<int>(cimg_image._height)-1 ) )
		return 0;

	//Linear interpolation
	double neighbor0 = cimg_image(xPixelValue, lowerPixelIndex );
	double neighbor1 = cimg_image(xPixelValue, higherPixelIndex );

	return neighbor0 + (neighbor1 - neighbor0) * (yCoordinateInPixel - lowerPixelIndex);
}

double Object2D::linear_atX(int yPixelValue, double xCoordinateInMM) const{
	/** Linear interpolation in X direction at xCoordinateInMM exactly at yPixelValue
	 *
	 */
	double xCoordinateInPixel = (objWidthHeight[0]/2 + xCoordinateInMM) / objPixSizes[0];
	int lowerPixelIndex = floor(xCoordinateInPixel);
	int higherPixelIndex = ceil(xCoordinateInPixel);
	if( ( lowerPixelIndex < 0 ) || ( higherPixelIndex > static_cast<int>(cimg_image._height)-1 ) )
		return 0;

	//Linear interpolation
	double neighbor0 = cimg_image(lowerPixelIndex, yPixelValue);
	double neighbor1 = cimg_image(higherPixelIndex, yPixelValue);

	return neighbor0 + (neighbor1 - neighbor0) * (xCoordinateInPixel - lowerPixelIndex);
}

double Object2D::getXValueAtPix(int pixIndex) const{
	/** Get the X value at the index pixIndex
	 *
	 */
	return xPixCentreCoords[pixIndex];
}

double Object2D::getYValueAtPix(int pixIndex) const{
	/** Get the Y value at the index pixIndex
	 *
	 */
	return yPixCentreCoords[pixIndex];
}

const Eigen::MatrixXd& Object2D::getDataAsEigenMatrixRef() const{
	/** Get the internal data structure as const reference
	 *
	 */
	return objData;
}
