#include <TomoLearnS/Gen1CT.hpp>
#include <iostream>
#include <cassert>
#include <cmath>

typedef std::vector<double> Row;
typedef std::vector<Row> Matrix;

Gen1CT::Gen1CT():detWidth{100},pixNum{100},object{nullptr}{
	};

Gen1CT::Gen1CT(double detWidth, int pixNum):detWidth{detWidth},pixNum{pixNum},object{nullptr} {
	pixPositions.resize(pixNum);
	double pixelDistance{detWidth/pixNum};
	for (int i=0; i<pixNum; i++){
		pixPositions[i]=-1*detWidth/2+(i+0.5)*pixelDistance;
	}
};

void Gen1CT::putObject(Object2D* sampleObject){
	object = sampleObject;
}

void Gen1CT::measure(const std::vector<double>& angles, int raysPerPixel){
	assert(object!=nullptr);

	const int numAngles = angles.size();

	Matrix sinogram(pixNum, Row(numAngles, 0)); //Matrix with (pixNum)x(numAngles) size, init with zeros

	cimg_library::CImg<uint16_t> sinoImage(pixNum, numAngles, 1, 1);

	double t, theta;
	for(int pixI=0; pixI<pixNum; ++pixI){
		t=pixPositions[pixI];
		for(int angI=0; angI<numAngles; ++angI){
			theta = std::fmod(angles[angI], M_PI);
			double sinTheta = sin(theta);
			double cosTheta = cos(theta);

			if( (theta > M_PI/4 ) && (theta < 3*M_PI/4)  ){
				for(int objectXIndex=0; objectXIndex < object -> getNumberOfPixels()[0]; ++objectXIndex){
					double objectYinMM = t*sinTheta+ (t*cosTheta - object->getXValueAtPix(objectXIndex))/sinTheta*cosTheta;
					double tmp = object->linear_atY(objectXIndex, objectYinMM) / std::abs(sinTheta)*object->getPixSizes()[0];
					sinogram[pixI][angI] += tmp;
				}
			} else{
				for(int objectYIndex=0; objectYIndex < object -> getNumberOfPixels()[1]; ++objectYIndex){
					double objectXinMM = t*cosTheta - (object->getYValueAtPix(objectYIndex)-t*sinTheta)/cosTheta*sinTheta;
					double tmp= object->linear_atX(objectYIndex, objectXinMM) / std::abs(cosTheta)*object->getPixSizes()[1];
					sinogram[pixI][angI] += tmp;
				}
			}

		}
	}

	double maxInt=0;
	for(int i=0; i<pixNum; ++i){
		for(int j=0; j<numAngles; ++j){
			if(sinogram[i][j] > maxInt) maxInt = sinogram[i][j];
		}
	}
	for(int i=0; i<pixNum; ++i){
		for(int j=0; j<numAngles; ++j){
			sinoImage(i,j) = static_cast<uint16_t>(sinogram[i][j]/maxInt*65536);
		}
	}

	cimg_library::CImgDisplay sinoWindow(sinoImage, "Sinogram");
	sinoWindow.wait();

}
