#include <Object2D.hpp>

#include <cstdint>

#ifdef Success       //Because otherwise Eigen not compile
  #undef Success
#endif

#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <thread>

/***
 * Constructor of Object2D which reads the data from a file on the disk
 * @param imageFilePath Path to the image file
 * @param objPixSizes X and Y size of a pixel
 * @param convertFromHUtoLA Convert data stored in HU to Linear Attenuation units
 * @param muWater //Linear attenuation coeff. of water [1/mm]
 */
Object2D::Object2D(const std::string& imageFilePath,
				   const std::array<double, 2>& objPixSizes,
				   bool convertFromHUtoLA,
				   double muWater):
															objPixSizes{objPixSizes} // @suppress("Symbol is not resolved")
				   {
	//Read the image file
    cimg_library::CImg<uint16_t> imageLoader(imageFilePath.c_str());
	//Copy data to Eigen::MatrixXd
	objData = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(imageLoader._width, imageLoader._height);
	for(uint i=0; i<imageLoader._width; ++i){
		for(uint j=0; j<imageLoader._height; ++j){
			objData(i,j) = imageLoader(i,j);
		}
	}

	//Convert Hounsfield (HU) to linear attenuation (LA) units
	if(convertFromHUtoLA){
		objData = (objData.array()-1000.0) * (muWater/1000) + muWater;   // HU -> LA transform
	}

	numberOfPixels = {static_cast<int>(objData.rows()), static_cast<int>(objData.cols())};
	xPixCentreCoords = std::vector<double>( static_cast<size_t>( objData.rows() ) );
	yPixCentreCoords = std::vector<double>( static_cast<size_t>( objData.cols() ) );

	objWidthHeightInMM = std::array<double, 2>{numberOfPixels[0] * objPixSizes[0], numberOfPixels[1] * objPixSizes[1]};

	std::cout << std:: endl << "Object in \"" << imageFilePath << "\" loaded" << std::endl <<
			"W*H [mm*mm]: " << objWidthHeightInMM[0] << " * "<< objWidthHeightInMM[1]<< std::endl;

	for(unsigned int i=0; i<imageLoader._width; ++i){
		xPixCentreCoords[i]=-1*objWidthHeightInMM[0]/2 + i*objPixSizes[0] + objPixSizes[0]/2;
	}
	for(unsigned int i=0; i<imageLoader._height; ++i){
		yPixCentreCoords[i]=   objWidthHeightInMM[1]/2 - i*objPixSizes[1] - objPixSizes[1]/2;
	}
}

/***
 * Constructor of Object2D which initialize the data with zeros
 * @param numberOfPixels Number of pixels of the data
 * @param objPixSizes Size os the pixels in X and Y directions
 */
Object2D::Object2D(const std::array<int, 2>& numberOfPixels,
		           const std::array<double, 2>& objPixSizes): numberOfPixels{numberOfPixels}, // @suppress("Symbol is not resolved")
		        		                                      objPixSizes{objPixSizes} { // @suppress("Symbol is not resolved")

	objData = Eigen::MatrixXd::Zero(numberOfPixels[0], numberOfPixels[1]);

	objWidthHeightInMM[0]=numberOfPixels[0]*objPixSizes[0];
    objWidthHeightInMM[1]=numberOfPixels[1]*objPixSizes[1];

    for(int i=0; i<numberOfPixels[0]; ++i){
    	xPixCentreCoords.push_back(-1*objWidthHeightInMM[0]/2 + i*objPixSizes[0] + objPixSizes[0]/2);
    }
    for(int i=0; i<numberOfPixels[1]; ++i){
    	yPixCentreCoords.push_back(objWidthHeightInMM[1]/2 - i*objPixSizes[1] - objPixSizes[1]/2);
    }
}

/***
 * Constructor which initializes Object2D with data from an Eigen::MatrixXd, detector width and vector of angles.
 * @param inData Data used for the initialization
 * @param detWidth Width of the detector array
 * @param angles Projection angles
 */
Object2D::Object2D(const Eigen::MatrixXd& inData,
				   double detWidth,
				   const Eigen::VectorXd& angles){
	objData=inData;
	numberOfPixels[0]=objData.rows();
	numberOfPixels[1]=objData.cols();

	objWidthHeightInMM[0]=detWidth;
	objPixSizes[0]=detWidth/numberOfPixels[0];

	objWidthHeightInMM[1]=numberOfPixels[1]*(angles[1]-angles[0]);
	objPixSizes[1]=angles[1]-angles[0];         //Suppose that angles are distributed uniformly

	xPixCentreCoords.resize(static_cast<size_t>( numberOfPixels[0] ) );
	for(int i=0; i<numberOfPixels[0]; i++){
		xPixCentreCoords[static_cast<size_t>(i)]=-1*objWidthHeightInMM[0]/2+(i+0.5)*objPixSizes[0];
	}
	yPixCentreCoords.resize(static_cast<size_t>(numberOfPixels[1]));
	for(int i=0; i<numberOfPixels[1]; i++){
		yPixCentreCoords[static_cast<size_t>(i)]= angles[i];
	}
}

/***
 * Constructor which initializes Object2D with data from an Eigen::MatrixXd and the pixel sizes
 * @param inData Data used for the initialization
 * @param objPixSizes Size of the pixels in X and Y directions
 */
Object2D::Object2D(const Eigen::MatrixXd& inData, const std::array<double, 2>& objPixSizes):objData{inData},
		                                                                                    objPixSizes{objPixSizes}{ // @suppress("Symbol is not resolved")
	numberOfPixels[0]=inData.rows();
	numberOfPixels[1]=inData.cols();

	objWidthHeightInMM[0]=objPixSizes[0]*numberOfPixels[0];
	objWidthHeightInMM[1]=objPixSizes[1]*numberOfPixels[1];

	xPixCentreCoords.resize(static_cast<size_t>( numberOfPixels[0] ) );
	for(int i=0; i<numberOfPixels[0]; i++){
		xPixCentreCoords[static_cast<size_t>(i)]=-1*objWidthHeightInMM[0]/2+(i+0.5)*objPixSizes[0];
	}
	yPixCentreCoords.resize(static_cast<size_t>( numberOfPixels[1] ) );
	for(int i=0; i<numberOfPixels[1]; i++){
		yPixCentreCoords[static_cast<size_t>(i)]=1*objWidthHeightInMM[1]/2-(i+0.5)*objPixSizes[1];
	}
}
/***
 * Destructor of Object2D class
 */
Object2D::~Object2D(){
	if(!cimg_window.is_closed()){
	    cimg_window.close();
	    cimg_window.flush();
	}
	if(displayThread.joinable())
		displayThread.join();
}

/**
 * Copy constructor of Object2D class which copies everything except the display which is bound to the original
 * @param objToCopy Object which is copied
 */
Object2D::Object2D(const Object2D& objToCopy){
	objData = objToCopy.objData;
	numberOfPixels = objToCopy.numberOfPixels;
	objPixSizes = objToCopy.objPixSizes;
	objWidthHeightInMM = objToCopy.objWidthHeightInMM;

	xPixCentreCoords = objToCopy.xPixCentreCoords;
	yPixCentreCoords = objToCopy.yPixCentreCoords;
}

/***
 * Copy assignment operator of Object2D which copies everything except the display which is bound to the original
 * @param objToCopy
 * @return Reference to the copied object
 */
Object2D& Object2D::operator=(const Object2D& objToCopy){
	objData = objToCopy.objData;
	numberOfPixels = objToCopy.numberOfPixels;
	objPixSizes = objToCopy.objPixSizes;
	objWidthHeightInMM = objToCopy.objWidthHeightInMM;

	xPixCentreCoords = objToCopy.xPixCentreCoords;
	yPixCentreCoords = objToCopy.yPixCentreCoords;

	return *this;
}

/***
 * Move constructor of Object2D moves everything except the display
 * @param objToMove Object to move
 */
Object2D::Object2D(Object2D&& objToMove):objData{std::move(objToMove.objData)},
                                         numberOfPixels{std::move(objToMove.numberOfPixels)}, // @suppress("Symbol is not resolved")
										 objPixSizes{std::move(objToMove.objPixSizes)}, // @suppress("Symbol is not resolved")
										 objWidthHeightInMM{std::move(objToMove.objWidthHeightInMM)}, // @suppress("Symbol is not resolved")
										 xPixCentreCoords{std::move(objToMove.xPixCentreCoords)}, // @suppress("Symbol is not resolved")
										 yPixCentreCoords{std::move(objToMove.yPixCentreCoords)}{ // @suppress("Symbol is not resolved")

}

/***
 * Move assignment operator of Object2D moves everything except the display
 * @param objToMove Object to move
 * @return Reference to the new object
*/
Object2D& Object2D::operator=(Object2D&& objToMove) noexcept{
	objData = std::move(objToMove.objData);
	numberOfPixels = std::move(objToMove.numberOfPixels);
	objPixSizes = std::move(objToMove.objPixSizes);
	objWidthHeightInMM = std::move(objToMove.objWidthHeightInMM);
	xPixCentreCoords = std::move(objToMove.xPixCentreCoords);
	yPixCentreCoords = std::move(objToMove.yPixCentreCoords);

	return *this;
}

/***
 * Add a scalar to the elements of matrix
 * @param addVal The value which is added to every element
 * @return Object2D with the new values
 */
Object2D Object2D::operator+(double addVal) const{
	Object2D result = *this;
	result.objData.array() += addVal;
	return result;
}

/***
 * Subtract a scalar from the elements of matrix
 * @param subVal The value which is subtracted from every element
 * @return Object2D with the new values
 */
Object2D Object2D::operator-(double subVal) const{
	return *this+(-1.0*subVal);
}

/***
 * Multiply every element of the matrix with the same scalar
 * @param coeff Scalar used for the multiplication
 * @return Object2D with the new values
 */
Object2D Object2D::operator*(double coeff) const{
	Object2D result = *this;
	result.objData = result.objData * coeff;
	return result;
}

/***
 * Set the value of a date point
 * @param i Row index
 * @param j Column index
 * @param setValue New value
 */
void Object2D::setData(int i, int j, double setValue){
	objData(i,j)=setValue;
}

/***
 * Display the date in a separate, interactive window.
 * @param label Label appearing in the window title
 */
void Object2D::display(const std::string& label){
	if(!(cimg_window.is_closed())){
		return;
	}
	if(displayThread.joinable()){
		displayThread.join();
	}

	std::stringstream ss;
	ss << " Label: " << label << " " << numberOfPixels[0] << " x " << numberOfPixels[1] << " points; "
			                         << std::fixed << std::setprecision(2)
		                             << objPixSizes[0] << " x " << objPixSizes[1] << " pixSize;"
						             << numberOfPixels[0]*objPixSizes[0] << " x "
									 << numberOfPixels[1]*objPixSizes[1] << " w x h;";
	std::string longTitle = ss.str();

	cimg_image.assign(static_cast<unsigned int>(numberOfPixels[0]), static_cast<unsigned int>(numberOfPixels[1]), 1, 1);

	for(unsigned int i=0; i<static_cast<unsigned int>(numberOfPixels[0]); ++i){
		for(unsigned int j=0; j<static_cast<unsigned int>(numberOfPixels[1]); ++j){
			//Without normalization
			cimg_image(i,j) = objData(i,j);
		}
	}

	char windowTitle[1024];
	if(longTitle.length() <=1024)
		strcpy(windowTitle, longTitle.c_str());
	else{
		longTitle.resize(1024);
		strcpy(windowTitle, longTitle.c_str());
	}

	displayThread = std::thread([this, windowTitle](){cimg_image._display(cimg_window,windowTitle,false,0,false,false); });
}

/***
 * Return the size of pixels in units of the coordinate system
 * @return A copy of the objPixSizes variable
 */
std::array<double,2> Object2D::getPixSizes() const{
	return objPixSizes;
}

/***
 * Returns the number of pixels
 * @return A copy of the numberOfPixels variable
 */
std::array<int, 2> Object2D::getNumberOfPixels() const{
	return numberOfPixels;
}

/***
 * Get the internal data structure as const reference
 * @return A const reference to the objData variable
 */
const Eigen::MatrixXd& Object2D::getDataAsEigenMatrixRef() const{
	return objData;
}

/**
 * Get the positions of pixel centers in X direction
 * @return Copy of the xPixCentreCoords variable
 */
std::vector<double> Object2D::getXPixCentreCoords() const{
	return xPixCentreCoords;
}

/**
 * Get the positions of pixel centers in Y direction
 * @return Copy of the yPixCentreCoords variable
 */
std::vector<double> Object2D::getYPixCentreCoords() const{
	return yPixCentreCoords;
}

/***
 * @brief Compares the data with the reference data. Returns the L2 norm.
 *
 * Performs bilinear interpolation on the compData data to match the coordinates
 * @param compData Data that is used for compare
 * @return The calculated L2 norm.
 */
double Object2D::compareNorm(Object2D compData) const{

	double L2norm=0;

	for(int xIdx=0; xIdx < numberOfPixels[0]; ++xIdx){
		for(int yIdx=0; yIdx < numberOfPixels[1]; ++yIdx){
			//Calculate the coordinates of the current pixel in the other image
			const double x = xPixCentreCoords[static_cast<size_t>(xIdx)];
			const double y = yPixCentreCoords[static_cast<size_t>(yIdx)];
			double xCoord = ( x + compData.objWidthHeightInMM[0]/2)/compData.objPixSizes[0]; // [pixel] in compData
			double yCoord = (compData.objWidthHeightInMM[1]/2 - y )/compData.objPixSizes[1]; // [pixel] in compData

			//Interpolate from the other image
			if( xCoord>=0.5 && xCoord <= numberOfPixels[0]-0.5 && yCoord >=0.5 && yCoord <= numberOfPixels[1]-0.5){
				double x1 = compData.xPixCentreCoords[floor(xCoord-0.5)];
				double x2 = compData.xPixCentreCoords[ceil(xCoord-0.5)];
				double y1 = compData.yPixCentreCoords[ceil(yCoord-0.5)];
				double y2 = compData.yPixCentreCoords[floor(yCoord-0.5)];

				double Q11=compData.objData(floor(xCoord), ceil(yCoord));
				double Q21=compData.objData(ceil(xCoord), ceil(yCoord));
				double Q12=compData.objData(floor(xCoord), floor(yCoord));
				double Q22=compData.objData(ceil(xCoord), floor(yCoord));

				double interpVal=1/compData.objPixSizes[0]/compData.objPixSizes[1]*
						                 (Q11*(x2-x)*(y2-y)+
										  Q21*(x-x1)*(y2-y)+
										  Q12*(x2-x)*(y-y1)+
										  Q22*(x-x1)*(y-y1));
				//Calculate the L2 norm
				L2norm+=std::pow(interpVal-objData(xIdx,yIdx),2);
			}
		}
	}

	//Normlalize the L2 norm
	L2norm /= numberOfPixels[0]*numberOfPixels[1];
	return L2norm;
}
