#pragma once
#include <string>
#include <Eigen/Dense>

class CTScan{
public:
	CTScan(std::string scanID, int pixNum, double detWidth, int numAngles);     //REGI
private:
	const std::string scanID;    //REGI
	const int pixNum;    //REGI
	const double detWidth;    //REGI
	const int numAngles;    //REGI
	Eigen::MatrixXd sinogram;    //REGI
};
