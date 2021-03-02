#pragma once
#include <TomoLearnS/Object2D.hpp>
#include <string>
#include <thread>
#include <Eigen/Dense>

#include <iostream>

class CTScan : public Object2D{
public:
	CTScan(std::string scanID, Eigen::MatrixXd sinogram, double detWidth, const Eigen::VectorXd& angles);

	const Eigen::VectorXd& getAnglesConstRef() const;

private:
	const std::string scanID;
	const Eigen::VectorXd angles;
};
