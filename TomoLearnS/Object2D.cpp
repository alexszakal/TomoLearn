#include <TomoLearnS/Object2D.hpp>

#include <cstdint>

#include <CImg.h>
#ifdef Success       //Because otherwise Eigen not compile
  #undef Success
#endif

#include <matplotlibcpp/matplotlibcpp_old.h>


#include <iostream>
#include <vector>
#include <array>

Object2D::Object2D(const std::string& imageFilePath,
				   const std::array<double, 2>& objPixSizes):
															objPixSizes{objPixSizes}
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

	/*//DEBUG show values along lower ellipses
	Eigen::VectorXd slice = image.row(821);
	matplotlibcpp::plot(std::vector<float> (&slice[0], slice.data()+slice.cols()*slice.rows()) );
	matplotlibcpp::show(False);
*/
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
    	yPixCentreCoords.push_back(objWidthHeightInMM[0]/2 - i*objPixSizes[1] - objPixSizes[1]/2);
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

	xPixCentreCoords.reserve(numberOfPixels[0]);
	for(int i=0; i<numberOfPixels[0]; i++){
		xPixCentreCoords[i]=-1*objWidthHeightInMM[0]/2+(i+0.5)*objPixSizes[0];
	}
	yPixCentreCoords.reserve(numberOfPixels[1]);
	for(int i=0; i<numberOfPixels[1]; i++){
			yPixCentreCoords[i]=-1*objWidthHeightInMM[1]/2+(i+0.5)*objPixSizes[1];
	}

}

/*Object2D::Object2D(Object2D&& objToMove){
	std::cout << "\nObject2D move constructor called!\n";
}*/

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
