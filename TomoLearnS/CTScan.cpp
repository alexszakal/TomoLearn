#include <TomoLearnS/CTScan.hpp>
#include  <TomoLearnS/Object2D.hpp>

    //REGI
CTScan::CTScan(std::string scanID, int pixNum, double detWidth, int numAngles):
					Object2D{std::array<int,2>{pixNum, numAngles},
	                std::array<double,2>{detWidth/pixNum, 1.0}},
					scanID{scanID}{

};

