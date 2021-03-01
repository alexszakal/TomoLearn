#pragma once
#include <TomoLearnS/Object2D.hpp>
#include <string>

class Reconst : public Object2D{
public:
	Reconst(std::string scanID, Eigen::MatrixXd recImage,
			const std::array<int, 2>& numberOfPixels, const std::array<double, 2>& objPixSizes);
	void display();
private:
	std::string scanID;
};



