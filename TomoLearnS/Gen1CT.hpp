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

	void measure(const std::string& label, const Eigen::VectorXd& angles, const std::string& scanLabel);
	void displayMeasurement(const std::string& label);

	void filteredBackProject(std::string sinogramID,
			                 const std::array<int,2>& numberOfRecPoints,
							 const std::array<double,2>& resolution );

	CTScan applyFilter(const std::string& sinogramID, Filter filter);
	void backProject(const CTScan& sinogram, const std::array<int,2>& numberOfRecPoints,
			                                 const std::array<double,2>& resolution);

private:
	const double detWidth;
	const int pixNum;
	std::vector<double> pixPositions;

	std::map<std::string, Phantom> phantoms;
	std::map<std::string, CTScan> scans;
	std::map<std::string, Reconst> reconsts;

};
