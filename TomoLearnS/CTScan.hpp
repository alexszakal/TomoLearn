#pragma once
#include  <TomoLearnS/Object2D.hpp>
#include <string>
#include <Eigen/Dense>

class CTScan : public Object2D{
public:
	CTScan(std::string scanID, int pixNum, double detWidth, int numAngles);     //REGI
private:
	const std::string scanID;    //REGI
	const int pixNum;    //REGI
	const double detWidth;    //REGI
	const int numAngles;    //REGI
	Eigen::MatrixXd sinogram;    //REGI
};
