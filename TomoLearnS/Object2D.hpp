#pragma once

#include <CImg.h>
#ifdef Success       //Because otherwise Eigen not compile
  #undef Success
#endif

#include <string>
#include <array>
#include <cstdint>
#include <vector>

typedef std::vector<double> Row ;
typedef  std::vector<Row> Matrix;

class Object2D{
public:
	Object2D():cimg_image{},cimg_window{},objPixSizes{1,1}{}  //Default constructor, initialize an empty image
	Object2D(const std::string& imageFilePath, const std::array<double, 2>& objPixSizes={0.1, 0.1});
	void display(std::string title);
	std::array<double, 2> getPixSizes();
	std::array<unsigned int, 2> getNumberOfPixels();
	double linear_atY(int xPixelValue, double yCoordinateInMM);
	double linear_atX(int yPixelValue, double xCoordinateInMM);
	double getXValueAtPix(int pixValue);
	double getYValueAtPix(int pixValue);
private:
	cimg_library::CImg<uint8_t> cimg_image;
	cimg_library::CImgDisplay cimg_window;
	const std::array<double, 2> objPixSizes;  //Size of a single pixel
	std::array<double, 2> objectSizeInMM;

	std::vector<double> xPixCentresInMM;  //
	std::vector<double> yPixCentresInMM;


};
