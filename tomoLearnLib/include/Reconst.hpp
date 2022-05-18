#pragma once
#include <Object2D.hpp>
#include <string>

class Reconst : public Object2D{
public:
	Reconst(): Object2D(), scanID("Empty"), convergenceCurve(std::vector<double>(0)) {
	}

	Reconst(std::string scanID,
			Eigen::MatrixXd recImage,
			const std::array<double, 2>& objPixSizes,
			std::vector<double> convergenceCurve=std::vector<double>(0)): Object2D(recImage, objPixSizes),
										   scanID(scanID),
										   convergenceCurve(convergenceCurve)
    {
    }

	std::vector<double> getConvergenceCurve();
private:
	std::string scanID;
	std::vector<double> convergenceCurve;
};



