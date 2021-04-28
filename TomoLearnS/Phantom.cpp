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

Phantom::Phantom( std::string label, const Object2D& dataPar):Object2D{dataPar},
		                                                     label{label}{
}

Phantom Phantom::operator*(double coeff) const {
	return Phantom{label, static_cast<Object2D>(*this) * coeff};
}

Phantom Phantom::operator+(double addVal) const {
	return Phantom{label, static_cast<Object2D>(*this) + addVal};
}

Phantom Phantom::operator-(double subVal) const {
	return *this + (-1.0*subVal);
}








