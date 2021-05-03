#pragma once

//#include <TomoLearnS/Object2D.hpp>
#include <TomoLearnS/CTScan.hpp>
#include <TomoLearnS/Phantom.hpp>
#include <TomoLearnS/Reconst.hpp>
#include <TomoLearnS/Filter.hpp>

#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

class Gen1CT{
public:
	Gen1CT();      //REGI
	Gen1CT(double detWidth, int pixNum);

	void addPhantom(const std::string& label,
			        const std::string& phantomImageSource,
					std::array<double, 2> pixSizes={0.1, 0.1});
	void displayPhantom(const std::string& label,
			            const std::string& title = "emptyTitle");

    void setI0(double newI0);

	void measure(const std::string& label, const Eigen::VectorXd& angles, const std::string& scanLabel);
	void displayMeasurement(const std::string& label);

	void filteredBackProject(std::string sinogramID,
			                 const std::array<int,2>& numberOfRecPoints,
							 const std::array<double,2>& resolution,
							 FilterType filterType,
							 double cutOffFreq,
							 const std::string& imageID);

	CTScan applyFilter(const std::string& sinogramID, Filter filter);
	Eigen::MatrixXd backProject(const CTScan& sinogram, const std::array<int,2>& numberOfRecPoints,
			                                 const std::array<double,2>& resolution);
	void displayReconstruction(const std::string& label);

	void compareRowPhantomAndReconst(int rowNum, const std::string& phantomID, const std::string& reconstID);
private:
	const double detWidth;
	const int pixNum;
	std::vector<double> pixPositions;

	double I0 = 3.0;   //Intensity of the tube

	double muWater = 0.02; //Linear attenuation coeff. of water [1/mm]!!!!

	std::map<std::string, Phantom> phantoms;
	std::map<std::string, CTScan> scans;
	std::map<std::string, Reconst> reconsts;

};
