#include <TomoLearnS/CTScan.hpp>
#include  <TomoLearnS/Object2D.hpp>

    //REGI
CTScan::CTScan(std::string scanID,  Eigen::MatrixXd sinogram,
		       int pixNum, double detWidth, const Eigen::VectorXd& angles):
					Object2D{scanID,
	                         std::array<int,2>{pixNum, angles.size()},
	                         std::array<double,2>{detWidth/pixNum, 1.0}
                            }
{

};

