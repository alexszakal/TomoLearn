#pragma once

#include <TomoLearnS/Object2D.hpp>
#include <vector>
#include <Eigen/Dense>

class Gen1CT{
public:
	Gen1CT();
	Gen1CT(double detWidth, int pixNum);
	void putObject(Object2D *object);
	void measure(const std::vector<double>& angles, int raysPerPixel=1);
	void filteredBackProjection();
	void displayMeasurement();

private:
	const double detWidth;
	const int pixNum;
	int numAngles;
	std::vector<double> pixPositions;
	Object2D* object;
	Eigen::MatrixXd sinogramE;
	cimg_library::CImg<uint16_t> sinoImage;
	cimg_library::CImgDisplay sinoWindow;


};
