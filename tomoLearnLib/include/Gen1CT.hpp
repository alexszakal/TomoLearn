#pragma once

//#include <TomoLearnS/Object2D.hpp>
#include <CTScan.hpp>
#include <Phantom.hpp>
#include <Reconst.hpp>
#include <Filter.hpp>

#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

#include <config.h>

enum class projectorType{pixelDriven,
                         Siddon,
                         rayDriven,
#if ENABLE_CUDA
						 rayDriven_GPU,
#endif
						 rayDrivenOptimized,
						 rayDrivenOptimizedFiniteBeam
                         };

enum class backprojectorType{pixelDriven,
#if ENABLE_CUDA
							 rayDriven_GPU,
#endif
                             rayDriven
                             };

enum class regularizerType{none,
						   quadratic,
						   Huber,
						   Gibbs
						   };

class Gen1CT{
public:
	Gen1CT();      //REGI
	Gen1CT(double detWidth, int pixNum);

	void addPhantom(const std::string& label,
			        const std::string& phantomImageSource,
					std::array<double, 2> pixSizes={0.1, 0.1},
					bool convertFromHUtoLA = false);

	void addPhantom(const Phantom& newPhantom);

	void displayPhantom(const std::string& label);

    void setI0(double newI0);

    void measure(const std::string& label, const Eigen::VectorXd& angles, const std::string& scanLabel, projectorType projector);

    Eigen::MatrixXd project(const Phantom& actualPhantom, const Eigen::VectorXd& angles, projectorType projector);

	void measure_withInterpolation(const std::string& label, const Eigen::VectorXd& angles, const std::string& scanLabel);

	Eigen::MatrixXd project_pixelDriven_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

	void measure_Siddon(const std::string& label, const Eigen::VectorXd& angles, const std::string& scanLabel);

	Eigen::MatrixXd project_Siddon_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

	Eigen::MatrixXd project_rayDriven_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

	Eigen::MatrixXd project_rayDrivenOptimized_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

	Eigen::MatrixXd project_rayDrivenOptimized_finiteBeam_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

#if ENABLE_CUDA
	Eigen::MatrixXd project_rayDriven_GPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);
#endif

	void displayMeasurement(const std::string& label);

	CTScan getMeasurement(const std::string& label);

	void filteredBackProject(std::string sinogramID,
			                 const std::array<int,2>& numberOfRecPoints,
							 const std::array<double,2>& resolution,
							 FilterType filterType,
							 double cutOffFreq,
							 backprojectorType backProjectAlgo,
							 const std::string& imageID);

	void MLEMReconst(std::string sinogramID,
				     const std::array<int,2>& numberOfRecPoints,
					 const std::array<double,2>& resolution,
					 projectorType projectAlgo,
					 backprojectorType backProjectAlgo,
					 const std::string& imageID,
					 int numberOfIterations);

	void SPSReconst(std::string sinogramID,
				     const std::array<int,2>& numberOfRecPoints,
					 const std::array<double,2>& resolution,
					 projectorType projectAlgo,
					 backprojectorType backProjectAlgo,
					 const std::string& imageID,
					 int numberOfIterations,
					 regularizerType regularizerFunction,
					 double beta=1000.0,
					 double delta=0.004);

	CTScan applyFilter(const std::string& sinogramID, Filter filter);

	Eigen::MatrixXd backProject(const CTScan& sinogram,
										const std::array<int,2>& numberOfRecPoints,
										const std::array<double,2>& resolution,
										backprojectorType backProjector);

	Eigen::MatrixXd backProject_pixelDriven_CPU(const CTScan& sinogram, const std::array<int,2>& numberOfRecPoints,
			                                 const std::array<double,2>& resolution);

	Eigen::MatrixXd backProject_rayDriven_CPU(const CTScan& sinogram, const std::array<int,2>& numberOfRecPoints,
				                                 const std::array<double,2>& resolution);

#if ENABLE_CUDA
	Eigen::MatrixXd backProject_rayDriven_GPU(const CTScan& sinogram, const std::array<int,2>& numberOfRecPoints,
					                                 const std::array<double,2>& resolution);
#endif

	void displayReconstruction(const std::string& label);

	void compareRowPhantomAndReconst(char direction, double position, const std::string& phantomID, const std::string& reconstID);

	void printPhantomParams(const std::string& phantomLabel);

#if ENABLE_CUDA
	void printGpuParameters();
#endif

private:
	const double detWidth;
	const int pixNum;
	std::vector<double> pixPositions; ///Positon of the pixel centers

	double I0 = 3.0;   //Intensity of the tube

	double muWater = 0.02; //Linear attenuation coeff. of water [1/mm]!!!!

	std::map<std::string, Phantom> phantoms;
	std::map<std::string, CTScan> scans;
	std::map<std::string, Reconst> reconsts;

};
