#pragma once
#include <Object2D.hpp>
#include <string>

class Reconst : public Object2D{
public:
	Reconst(): Object2D(), scanID("Empty") {
	}

	Reconst(std::string scanID,
			Eigen::MatrixXd recImage,
			const std::array<double, 2>& objPixSizes): Object2D(recImage,
																objPixSizes),
													   scanID(scanID)
    {
    }
private:
	std::string scanID;
};



