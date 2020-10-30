#pragma once
#include <string>
#include <Eigen/Dense>

class CTScan{
public:
	CTScan(std::string scanID, int pixNum, double detWidth, int numAngles);
private:
	const std::string scanID;
	const int pixNum;
	const double detWidth;
	const int numAngles;
	Eigen::MatrixXd sinogram;
};
