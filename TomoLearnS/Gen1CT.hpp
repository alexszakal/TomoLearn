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
	void backProject(const std::vector<int>& numberOfRecPoints, const std::vector<double>& resolution);
	void FBP(std::vector<int> numberOfRecPoints, std::vector<double> resolution);
	void FBP_bandlimited(std::vector<int> numberOfRecPoints, std::vector<double> resolution);
	const Eigen::MatrixXd& getSinogram() const;
	const std::vector<double>& getPixPositions() const;

private:
	const double detWidth;
	const int pixNum;
	Object2D* object;
	int numAngles;
	std::vector<double> angs;
	std::vector<double> pixPositions;

	Eigen::MatrixXd sinogram;
	cimg_library::CImg<uint16_t> sinoImage;
	cimg_library::CImgDisplay sinoWindow;

	Eigen::MatrixXd backprojection;
	cimg_library::CImg<uint16_t> BPImage;
	cimg_library::CImgDisplay BPWindow;

	Eigen::MatrixXd filteredBP;
	cimg_library::CImg<uint16_t> FBPImage;
	cimg_library::CImgDisplay FBPWindow;

};
