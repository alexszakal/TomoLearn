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
						 rayDrivenOptimized
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
/**
 * @brief Gen1CT class implements the projectors, backprojectors, reconstruction algorithms and data storage.
 */
class Gen1CT{
public:
	Gen1CT();
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

	Eigen::MatrixXd project_pixelDriven_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

	Eigen::MatrixXd project_Siddon_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

	Eigen::MatrixXd project_rayDriven_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

	Eigen::MatrixXd project_rayDrivenOptimized_CPU(const Phantom& actualPhantom, const Eigen::VectorXd& angles);

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
							 const std::string& imageID,
							 std::string referenceImage="");

	void MLEMReconst(std::string sinogramID,
				     const std::array<int,2>& numberOfRecPoints,
					 const std::array<double,2>& resolution,
					 projectorType projectAlgo,
					 backprojectorType backProjectAlgo,
					 const std::string& imageID,
					 int numberOfIterations,
					 std::string referenceImage="");

	void SPSReconst(std::string sinogramID,
				     const std::array<int,2>& numberOfRecPoints,
					 const std::array<double,2>& resolution,
					 projectorType projectAlgo,
					 backprojectorType backProjectAlgo,
					 const std::string& imageID,
					 int numberOfIterations,
					 regularizerType regularizerFunction,
					 double beta=2000.0,
					 double delta=0.004,
					 std::string referenceImage="");

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

	double compareReconToPhantom(std::string reconLabel, std::string phantomLabel) const;

	std::vector<double> getConvergenceCurve(std::string label) const;

private:
	const double detWidth;  /** Detector width [mm] */
	const int pixNum; /** Number of detector pixels [1] */
	std::vector<double> pixPositions; /** Positon of the pixel centers */

	double I0 = 3.0;   /** Intensity of the tube */

	double muWater = 0.02; /** Linear attenuation coeff. of water [1/mm] */

	std::map<std::string, Phantom> phantoms;  /** Library storing the Phantoms*/
	std::map<std::string, CTScan> scans;      /** Library storing the Scans*/
	std::map<std::string, Reconst> reconsts;  /** Library storing the Reconstructions*/
};
