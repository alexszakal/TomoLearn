#include <TomoLearnS/CTScan.hpp>
#include  <TomoLearnS/Object2D.hpp>

#include <iostream>

CTScan::CTScan(std::string scanID,
		       Eigen::MatrixXd sinogram,
			   double detWidth,
			   const Eigen::VectorXd& angles):
				        	                 Object2D(sinogram, detWidth, angles),
											 scanID{scanID},
											 angles{angles}
{

}

const Eigen::VectorXd& CTScan::getAnglesConstRef() const{
	return angles;
}

