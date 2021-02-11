#pragma once
#include  <TomoLearnS/Object2D.hpp>
#include <string>
#include <Eigen/Dense>

class CTScan : public Object2D{
public:
	CTScan(std::string scanID, Eigen::MatrixXd sinogram, int pixNum, double detWidth, const Eigen::VectorXd& angles);
private:
	const std::string scanID;
};
