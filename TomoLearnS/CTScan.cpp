#include <TomoLearnS/CTScan.hpp>

    //REGI
CTScan::CTScan(std::string scanID, int pixNum, double detWidth, int numAngles):
					scanID{scanID}, pixNum{pixNum}, detWidth{detWidth},numAngles{numAngles}{
    sinogram =Eigen::MatrixXd::Zero(pixNum, numAngles);
};

