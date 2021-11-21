#pragma once
#include <Object2D.hpp>
#include <string>
#include <thread>
#include <Eigen/Dense>

#include <iostream>

class CTScan : public Object2D{
public:
	CTScan();
	CTScan(std::string scanID, Eigen::MatrixXd sinogram, double detWidth, const Eigen::VectorXd& angles);
	CTScan(const std::string& scanID, const Eigen::VectorXd& angles, const Object2D& dataPar);
	CTScan(std::string scanID,
			double detWidth,
			int numDetPixels,
			const Eigen::VectorXd& angles,
			const std::vector<double>& rhos,
			const std::vector<double>& alphas,
			const std::vector< std::array<double,2> >& centers,
			const std::vector< std::array<double,2> >& axes );

	const Eigen::VectorXd& getAnglesConstRef() const;
	double getDetWidth() const;

	friend CTScan operator/(const CTScan& lhs, const CTScan& rhs);
	friend CTScan operator+(const CTScan& lhs, double rhs);

private:
	std::string scanID;
	Eigen::VectorXd angles;
};
