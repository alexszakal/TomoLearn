#pragma once
#include  <TomoLearnS/Object2D.hpp>
#include <string>
#include <thread>
#include <Eigen/Dense>

class CTScan : public Object2D{
public:
	CTScan(std::string scanID, Eigen::MatrixXd sinogram, double detWidth, const Eigen::VectorXd& angles);
	void display();
private:
	const std::string scanID;
	const Eigen::VectorXd angles;

	cimg_library::CImg<uint16_t> cimg_image;
	cimg_library::CImgDisplay cimg_window;
	std::thread displayThread;
};
