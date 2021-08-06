#pragma once
#include <TomoLearnS/Object2D.hpp>
#include <string>
#include <thread>
#include <Eigen/Dense>

#include <iostream>

class CTScan : public Object2D{
public:
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



private:
	const std::string scanID;
	const Eigen::VectorXd angles;
};
