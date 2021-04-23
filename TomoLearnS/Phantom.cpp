#include <TomoLearnS/Phantom.hpp>
#include <iostream>

//Phantom::Phantom():Object2D(){
//
//}

Phantom::Phantom(const std::string& label,
		         const std::string& imageFilePath,
				 const std::array<double, 2>& objPixSizes):
				                                          Object2D{imageFilePath, objPixSizes},
														  label{label}{
}

Phantom::Phantom(const std::string& label,
		         const Eigen::MatrixXd inData,
				 const std::array<double, 2> objPixSizes):
						                                  Object2D{inData, objPixSizes},
														  label{label}{
}

Phantom Phantom::operator*(double coeff) const {
	return Phantom{label, this->getDataAsEigenMatrixRef()*coeff, this->getPixSizes()};
}

Phantom Phantom::operator+(double addVal) const {
	return Phantom{label, this->getDataAsEigenMatrixRef().array() + addVal, this->getPixSizes()};
}






