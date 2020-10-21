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
		for(int angI=0; angI<numAngles; ++angI){
			t=pixPositions[pixI];
			theta = angles[angI];
			double sinTheta = sin(theta);
			double cosTheta = cos(theta);
			std::cout << "t=" << t << " theta=" << theta/M_PI*180;
			//Go through the object with X values and interpolate in Y
			for(int objectXIndex=0; objectXIndex < object -> getNumberOfPixels()[0]; ++objectXIndex){
				double objectYinMM = t*sinTheta+ (t*cosTheta - object->getXValueAtPix(objectXIndex))/sinTheta*cosTheta;
				std::cout << " " << objectYinMM;
				double tmp = object->linear_atY(objectXIndex, objectYinMM) /std::abs(sinTheta)*object->getPixSizes()[0];
				if( object->linear_atY(objectXIndex, objectYinMM) > 0)
					std::cout << "t: " << t << " Theta: " << theta/M_PI*180 << " objectXIndex: " << objectXIndex << std::endl;
				sinogram[pixI][ angI] += tmp;
			}
			std::cout << std::endl;
		}
	}

/*
	for(int i; i<numAngles; ++i){
		std::cout << sinogram[150][i] <<"  ";
		if (i%10 ==0) std::cout << "\n";
	}
*/

	cimg_library::CImgDisplay sinoWindow(sinoImage, "Sinogram");
	sinoWindow.wait();


}
