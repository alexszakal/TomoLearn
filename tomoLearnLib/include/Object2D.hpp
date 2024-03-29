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

/***
 * @brief Base class to hold 2D maps of data with axes information and display capabilities.
 */
class Object2D{
public:
	//Default constructor, initialize an empty image
	Object2D():numberOfPixels{0,0},objPixSizes{1,1},objWidthHeightInMM{0,0}{}
	//Constructor for init with data from file
	Object2D(const std::string& imageFilePath, const std::array<double, 2>& objPixSizes={0.1, 0.1},
			 bool convertFromHUtoLA = false, double muWater=0.02);
	//Constructor of a zero-initialized Object2D
	Object2D(const std::array<int, 2>& numberOfPixels, const std::array<double, 2>& objPixSizes);
	//Constructor for initializing with Eigen::Matrix for Sinogram
	Object2D(const Eigen::MatrixXd& inData, double detWidth, const Eigen::VectorXd& angles);
	//Constructor for initialization with Eigen::Matrix for Phantom
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

	//Set data
	void setData(int i, int j, double setValue);

	//Getter functions
	std::array<double, 2> getPixSizes() const;
	std::array<int, 2> getNumberOfPixels() const;
	std::vector<double> getXPixCentreCoords() const;
	std::vector<double> getYPixCentreCoords() const;
	const Eigen::MatrixXd& getDataAsEigenMatrixRef() const;

	//Functions for measurement
	double linear_atY(int xPixelValue, double yCoordinateInMM) const;
	double linear_atX(int yPixelValue, double xCoordinateInMM) const;
	double getXValueAtPix(int pixValue) const ;
	double getYValueAtPix(int pixValue) const;

	void display(const std::string& label="");

	double compareNorm(Object2D compData) const;

private:
	//Raw Data
	//1st index (row number) -> X direction
	//2nd index (col number) -> Y direction
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> objData; /** The that is stored */

	//Axes parameters
	std::array<int, 2> numberOfPixels;  /** Number of data pixels in X and Y directions */
	std::array<double, 2> objPixSizes;   /** Sizes of a single pixel */
	std::array<double, 2> objWidthHeightInMM; /** Total size of the data in X and Y directions */

	std::vector<double> xPixCentreCoords; /**Coordinates of pixel centers in X direction */
	std::vector<double> yPixCentreCoords; /**Coordinates of pixel centers in Y direction */

	//Display
	cimg_library::CImg<double> cimg_image;  /** cimg object used for the display */
	cimg_library::CImgDisplay cimg_window;  /** cimg_window used for the display */
	std::thread displayThread;              /** Separate thread where the diplay function runs */

};

/***
 * Linear interpolation in Y direction at yCoordinateInMM exactly at xPixelValue
 * @param xPixelValue The number of X pixel where the interpolation performed
 * @param yCoordinateInMM Y coordinate in [mm] where the interpolated value is requested
 * @return Interpolated value
 */
inline double Object2D::linear_atY(int xPixelValue, double yCoordinateInMM) const{

	double yCoordinateInPixel = ( (objWidthHeightInMM[1]/2 - objPixSizes[1]/2) - yCoordinateInMM) / objPixSizes[1]; //Hany pixelre vagyunk a 0. pixeltol

	if( ( yCoordinateInPixel >= 0.0 ) && ( yCoordinateInPixel < numberOfPixels[1]-1.0 ) ){
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

/***
 * Linear interpolation in X direction at xCoordinateInMM exactly at yPixelValue
 * @param yPixelValue The number of Y pixel where the interpolation performed
 * @param xCoordinateInMM X coordinate in [mm] where the interpolated value is requested
 * @return Interpolated value
 */
inline double Object2D::linear_atX(int yPixelValue, double xCoordinateInMM) const{

	double xCoordinateInPixel = (   (objWidthHeightInMM[0]/2 - objPixSizes[0]/2) + xCoordinateInMM  ) / objPixSizes[0];

	if( ( xCoordinateInPixel > 0.0 ) && ( xCoordinateInPixel < numberOfPixels[0]-1.0 ) ){
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

/***
 * Get the X value at the index pixIndex
 * @param pixIndex The index where the X-value is returned
 * @return The X-value at the pixIndex position
 */
inline double Object2D::getXValueAtPix(int pixIndex) const{
	return xPixCentreCoords[static_cast<size_t>(pixIndex)];
}

/***
 * Get the Y value at the index pixIndex
 * @param pixIndex The index where the Y-value is returned
 * @return The Y-value at the pixIndex position
 */
inline double Object2D::getYValueAtPix(int pixIndex) const{
	return yPixCentreCoords[static_cast<size_t>(pixIndex)];
}


