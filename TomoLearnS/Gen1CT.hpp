#pragma once

//#include <TomoLearnS/Object2D.hpp>
#include <TomoLearnS/CTScan.hpp>
#include <TomoLearnS/Phantom.hpp>
#include <TomoLearnS/Reconst.hpp>

#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

class Gen1CT{
public:
	Gen1CT();      //REGI
	Gen1CT(double detWidth, int pixNum);    //REGI

	void addPhantom(const std::string& label, const std::string& phantomImageSource, std::array<double, 2> pixSizes={0.1, 0.1});
	void displayPhantom(const std::string& label, const std::string& title = "emptyTitle");

	void measure(const std::string& label, const Eigen::VectorXd& angles, const std::string& scanLabel);
	void displayMeasurement(const std::string& label);

	void filteredBackProjection();    //REGI
	void backProject(const std::vector<int>& numberOfRecPoints, const std::vector<double>& resolution);    //REGI
	void FBP(std::vector<int> numberOfRecPoints, std::vector<double> resolution);    //REGI

private:
	const double detWidth;
	const int pixNum;
	std::vector<double> pixPositions;

	std::map<std::string, Phantom> phantoms;
	std::map<std::string, CTScan> scans;
	std::map<std::string, Reconst> reconsts;

	//Eigen::MatrixXd sinogram;    //REGI
	//cimg_library::CImg<uint16_t> sinoImage;    //REGI
	//cimg_library::CImgDisplay sinoWindow;    //REGI

	//Eigen::MatrixXd backprojection;    //REGI
	//cimg_library::CImg<uint16_t> BPImage;    //REGI
	//cimg_library::CImgDisplay BPWindow;    //REGI

	//Eigen::MatrixXd filteredBP;    //REGI
	//cimg_library::CImg<uint16_t> FBPImage;    //REGI
	//cimg_library::CImgDisplay FBPWindow;    //REGI

};
