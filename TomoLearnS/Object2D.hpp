#pragma once

#include <CImg.h>
#ifdef Success       //Because otherwise Eigen not compile
  #undef Success
#endif
#include <Eigen/Dense>
#include <string>
#include <array>
#include <cstdint>
#include <vector>

class Object2D{
public:
	Object2D():objPixSizes{1,1},objWidthHeightInMM{0,0},numberOfPixels{0,0},cimg_image{},cimg_window{}{}  //Default constructor, initialize an empty image   //REGI
	Object2D(const std::string& label, const std::string& imageFilePath, const std::array<double, 2>& objPixSizes={0.1, 0.1});
	Object2D(const std::array<int, 2>& numberOfPixels, const std::array<double, 2>& objPixSizes);

	void display(const std::string& title="", bool isInterractive=true);
	std::array<double, 2> getPixSizes() const;
	std::array<int, 2> getNumberOfPixels() const;
	double linear_atY(int xPixelValue, double yCoordinateInMM) const;
	double linear_atX(int yPixelValue, double xCoordinateInMM) const;
	double getXValueAtPix(int pixValue) const ;
	double getYValueAtPix(int pixValue) const;
	const Eigen::MatrixXd& getDataAsEigenMatrixRef() const;
private:
	//Raw Data
	//1st index (row number) -> X direction
	//2nd index (col number) -> Y direction
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> objData;
	//Axes parameters
	const std::array<double, 2> objPixSizes;   //Size of a single pixel
	std::array<double, 2> objWidthHeightInMM;
	std::array<int, 2> numberOfPixels;

	const std::string label;

	std::vector<double> xPixCentreCoords;
	std::vector<double> yPixCentreCoords;
	//Display
	cimg_library::CImg<uint16_t> cimg_image;
	cimg_library::CImgDisplay cimg_window;
};
