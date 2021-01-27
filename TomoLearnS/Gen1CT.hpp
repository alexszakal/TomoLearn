#pragma once

#include <TomoLearnS/Object2D.hpp>
#include <TomoLearnS/CTScan.hpp>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

class Gen1CT{
public:
	Gen1CT();      //REGI
	Gen1CT(double detWidth, int pixNum);    //REGI

	void addPhantom(const std::string& label, const std::string& phantomImageSource, std::array<double, 2> pixSizes={0.1, 0.1});
	void displayPhantom(const std::string& label);

	void measure(const std::string& label, const Eigen::VectorXd& angles, const std::string& scanLabel);
	void displayMeasurement();    //REGI

	void filteredBackProjection();    //REGI
	void backProject(const std::vector<int>& numberOfRecPoints, const std::vector<double>& resolution);    //REGI
	void FBP(std::vector<int> numberOfRecPoints, std::vector<double> resolution);    //REGI
	void FBP_bandlimited(std::vector<int> numberOfRecPoints, std::vector<double> resolution);    //REGI
	const Eigen::MatrixXd& getSinogram() const;    //REGI
	const std::vector<double>& getPixPositions() const;    //REGI

private:
	const double detWidth;    //REGI
	const int pixNum;    //REGI
	Object2D* object;    //REGI
	int numAngles;    //REGI

	std::map<std::string, Object2D> phantoms;
	std::map<std::string, CTScan> scans;

	Eigen::VectorXd angs;    //REGI
	std::vector<double> pixPositions;    //REGI

	Eigen::MatrixXd sinogram;    //REGI
	cimg_library::CImg<uint16_t> sinoImage;    //REGI
	cimg_library::CImgDisplay sinoWindow;    //REGI

	Eigen::MatrixXd backprojection;    //REGI
	cimg_library::CImg<uint16_t> BPImage;    //REGI
	cimg_library::CImgDisplay BPWindow;    //REGI

	Eigen::MatrixXd filteredBP;    //REGI
	cimg_library::CImg<uint16_t> FBPImage;    //REGI
	cimg_library::CImgDisplay FBPWindow;    //REGI

};
