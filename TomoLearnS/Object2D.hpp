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
#include <thread>

class Object2D{
public:
	//Default constructor, initialize an empty image
	Object2D():objPixSizes{1,1},objWidthHeightInMM{0,0},numberOfPixels{0,0}{}
	//Constructor for init with data from file
	Object2D(const std::string& imageFilePath, const std::array<double, 2>& objPixSizes={0.1, 0.1});
	//Constructor of a zero-initialized Object2D
	Object2D(const std::array<int, 2>& numberOfPixels, const std::array<double, 2>& objPixSizes);
	//Constructor for initializing with Eigen::Matrix
	Object2D(const Eigen::MatrixXd& inData, double detWidth, const Eigen::VectorXd& angles);
	//Constructor for initialization with Eigen::Matrix
	Object2D(const Eigen::MatrixXd& inData, const std::array<double, 2>& objPixSizes={0.1, 0.1});

	~Object2D();
	//Rule of 5 because custom destructor defined.
	Object2D(const Object2D& objToCopy); //Copy constructor
	Object2D& operator=(const Object2D& objToCopy);		//Copy assignment
	Object2D(Object2D&& objToMove); //Move constructor
	Object2D& operator=(Object2D&& objToMove) noexcept;	//Move assignment

	//Arithmetic operators
	Object2D operator+(double addVal) const;
	Object2D operator*(double coeff) const;
	Object2D operator-(double subVal) const;


	std::array<double, 2> getPixSizes() const;
	std::array<int, 2> getNumberOfPixels() const;

	const Eigen::MatrixXd& getDataAsEigenMatrixRef() const;

	//Functions for measurement
	double linear_atY(int xPixelValue, double yCoordinateInMM) const;
	double linear_atX(int yPixelValue, double xCoordinateInMM) const;
	double getXValueAtPix(int pixValue) const ;
	double getYValueAtPix(int pixValue) const;

	void display(const std::string& label);

private:
	//Raw Data
	//1st index (row number) -> X direction
	//2nd index (col number) -> Y direction
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> objData;

	//Axes parameters
	std::array<int, 2> numberOfPixels;
	std::array<double, 2> objPixSizes;   //Size of a single pixel
	std::array<double, 2> objWidthHeightInMM;

	std::vector<double> xPixCentreCoords;
	std::vector<double> yPixCentreCoords;

	//Display
	cimg_library::CImg<double> cimg_image;
	cimg_library::CImgDisplay cimg_window;
	std::thread displayThread;

};

inline double Object2D::linear_atY(int xPixelValue, double yCoordinateInMM) const{
	/** Linear interpolation in X direction at yCoordinateInMM exactly at xPixelValue
	 *
	 */
	double yCoordinateInPixel = (objWidthHeightInMM[1]/2 - yCoordinateInMM) / objPixSizes[1];

	if( ( yCoordinateInPixel > 0 ) && ( yCoordinateInPixel < numberOfPixels[1]-1 ) ){
		int lowerPixelIndex = static_cast<int>(yCoordinateInPixel); //Casting is MUCH faster than floor and res. is same because yCoordinateInPixel > 0

		//Linear interpolation
		double neighbor0 = objData(xPixelValue, lowerPixelIndex );     //Data access limits the throughput
		double neighbor1 = objData(xPixelValue, lowerPixelIndex+1 );   //Data access limits the throughput

		return neighbor0 + (neighbor1 - neighbor0) * (yCoordinateInPixel - lowerPixelIndex);
	}
	else{
		return 0;
	}
}

inline double Object2D::linear_atX(int yPixelValue, double xCoordinateInMM) const{
	/** Linear interpolation in X direction at xCoordinateInMM exactly at yPixelValue
	 *
	 */
	double xCoordinateInPixel = (objWidthHeightInMM[0]/2 + xCoordinateInMM) / objPixSizes[0];

	if( ( xCoordinateInPixel > 0 ) && ( xCoordinateInPixel < numberOfPixels[0]-1 ) ){
		int lowerPixelIndex = static_cast<int>(xCoordinateInPixel); //Casting is MUCH faster than floor and res. is same because xCoordinateInPixel > 0

		//Linear interpolation
		double neighbor0 = objData(lowerPixelIndex, yPixelValue);   //Data access limits the throughput
		double neighbor1 = objData(lowerPixelIndex+1, yPixelValue); //Data access limits the throughput

		return neighbor0 + (neighbor1 - neighbor0) * (xCoordinateInPixel - lowerPixelIndex);
	}
	else{
		return 0;
	}
}

inline double Object2D::getXValueAtPix(int pixIndex) const{
	/** Get the X value at the index pixIndex
	 *
	 */
	return xPixCentreCoords[pixIndex];
}

inline double Object2D::getYValueAtPix(int pixIndex) const{
	/** Get the Y value at the index pixIndex
	 *
	 */
	return yPixCentreCoords[pixIndex];
}


