#include <Object2D.hpp>

#include <cstdint>

#ifdef Success       //Because otherwise Eigen not compile
  #undef Success
#endif

#include <matplotlibcpp.h>

#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <thread>

Object2D::Object2D(const std::string& imageFilePath,
				   const std::array<double, 2>& objPixSizes,
				   bool convertFromHUtoLA):
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
		double muWater = 0.02; //Linear attenuation coeff. of water [1/mm]!!!!
		objData = (objData.array()-1000.0) * (muWater/1000) + muWater;   // HU -> LA transform
	}

	numberOfPixels = {static_cast<int>(objData.rows()), static_cast<int>(objData.cols())};
	xPixCentreCoords = std::vector<double>( objData.rows()  );
	yPixCentreCoords = std::vector<double>( objData.cols() );

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

/* Initialize an Object2D with zeros*/
Object2D::Object2D(const std::array<int, 2>& numberOfPixels,
		           const std::array<double, 2>& objPixSizes): objPixSizes{objPixSizes}, // @suppress("Symbol is not resolved")
															  numberOfPixels{numberOfPixels} { // @suppress("Symbol is not resolved")

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

//Initialize object with data in Eigen::MatrixXd
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

	xPixCentreCoords.resize(numberOfPixels[0]);
	for(int i=0; i<numberOfPixels[0]; i++){
		xPixCentreCoords[i]=-1*objWidthHeightInMM[0]/2+(i+0.5)*objPixSizes[0];
	}
	yPixCentreCoords.resize(numberOfPixels[1]);
	for(int i=0; i<numberOfPixels[1]; i++){
		yPixCentreCoords[i]= angles[i];
	}
}

Object2D::Object2D(const Eigen::MatrixXd& inData, const std::array<double, 2>& objPixSizes):objData{inData},
		                                                                                    objPixSizes{objPixSizes}{ // @suppress("Symbol is not resolved")
	numberOfPixels[0]=inData.rows();
	numberOfPixels[1]=inData.cols();

	objWidthHeightInMM[0]=objPixSizes[0]*numberOfPixels[0];
	objWidthHeightInMM[1]=objPixSizes[1]*numberOfPixels[1];

	xPixCentreCoords.resize(numberOfPixels[0]);
	for(int i=0; i<numberOfPixels[0]; i++){
		xPixCentreCoords[i]=-1*objWidthHeightInMM[0]/2+(i+0.5)*objPixSizes[0];
	}
	yPixCentreCoords.resize(numberOfPixels[1]);
	for(int i=0; i<numberOfPixels[1]; i++){
		yPixCentreCoords[i]=1*objWidthHeightInMM[1]/2-(i+0.5)*objPixSizes[1];
	}
}

Object2D::~Object2D(){
	//destructor
	if(!cimg_window.is_closed()){
	    cimg_window.close();
	    cimg_window.flush();
	}
	if(displayThread.joinable())
		displayThread.join();
}

Object2D::Object2D(const Object2D& objToCopy){
	//Copy constructor copies everything except the display which is bound to the original
	objData = objToCopy.objData;
	numberOfPixels = objToCopy.numberOfPixels;
	objPixSizes = objToCopy.objPixSizes;   //Size of a single pixel
	objWidthHeightInMM = objToCopy.objWidthHeightInMM;

	xPixCentreCoords = objToCopy.xPixCentreCoords;
	yPixCentreCoords = objToCopy.yPixCentreCoords;
}

Object2D& Object2D::operator=(const Object2D& objToCopy){
	//Copy assignment copies everything except the display which is bound to the original
	objData = objToCopy.objData;
	numberOfPixels = objToCopy.numberOfPixels;
	objPixSizes = objToCopy.objPixSizes;   //Size of a single pixel
	objWidthHeightInMM = objToCopy.objWidthHeightInMM;

	xPixCentreCoords = objToCopy.xPixCentreCoords;
	yPixCentreCoords = objToCopy.yPixCentreCoords;

	return *this;
}

Object2D::Object2D(Object2D&& objToMove):objData{std::move(objToMove.objData)},
                                         numberOfPixels{std::move(objToMove.numberOfPixels)}, // @suppress("Symbol is not resolved")
										 objPixSizes{std::move(objToMove.objPixSizes)}, // @suppress("Symbol is not resolved")
										 objWidthHeightInMM{std::move(objToMove.objWidthHeightInMM)}, // @suppress("Symbol is not resolved")
										 xPixCentreCoords{std::move(objToMove.xPixCentreCoords)}, // @suppress("Symbol is not resolved")
										 yPixCentreCoords{std::move(objToMove.yPixCentreCoords)}{ // @suppress("Symbol is not resolved")
	//Move constructor
}

Object2D& Object2D::operator=(Object2D&& objToMove) noexcept{
	//Move assignment
	objData = std::move(objToMove.objData);
	numberOfPixels = std::move(objToMove.numberOfPixels);
	objPixSizes = std::move(objToMove.objPixSizes);
	objWidthHeightInMM = std::move(objToMove.objWidthHeightInMM);
	xPixCentreCoords = std::move(objToMove.xPixCentreCoords);
	yPixCentreCoords = std::move(objToMove.yPixCentreCoords);

	return *this;
}

Object2D Object2D::operator+(double addVal) const{
	Object2D result = *this;
	result.objData.array() += addVal;
	return result;
}

Object2D Object2D::operator-(double subVal) const{
	return *this+(-1.0*subVal);
}

Object2D Object2D::operator*(double coeff) const{
	Object2D result = *this;
	result.objData = result.objData * coeff;
	return result;
}

void Object2D::setData(int i, int j, double setValue){
	objData(i,j)=setValue;
}

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



const Eigen::MatrixXd& Object2D::getDataAsEigenMatrixRef() const{
	/** Get the internal data structure as const reference
	 *
	 */
	return objData;
}

std::vector<double> Object2D::getXPixCentreCoords() const{
	return xPixCentreCoords;
}

std::vector<double> Object2D::getYPixCentreCoords() const{
	return yPixCentreCoords;
}

/***
 * Compares the data with the reference data. Returns the L2 norm.
 * Performs bilinear interpolation on the compData data to match the coordinates
 * @param compData Data that is used for compare
 * @return The calculated L2 norm.
 */
double Object2D::compareNorm(Object2D compData) const{

	double L2norm=0;

	for(int xIdx=0; xIdx < numberOfPixels[0]; ++xIdx){
		for(int yIdx=0; yIdx < numberOfPixels[1]; ++yIdx){
			//Calculate the coordinates of the current pixel in the other image
			const double x = xPixCentreCoords[xIdx];
			const double y = yPixCentreCoords[yIdx];
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
