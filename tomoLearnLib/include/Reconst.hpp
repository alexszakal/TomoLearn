#pragma once
#include <Object2D.hpp>
#include <string>

/***
 * @brief Class to store the reconstructions
 */
class Reconst : public Object2D{
public:
	/***
	 * Constructor of an empty Reconst object
	 */
	Reconst(): Object2D(), scanID("Empty"), convergenceCurve(std::vector<double>(0)) {
	}

	/***
	 * Construct a REconst object from a matrix of the reconstructed data
	 * @param scanID Identifier label of the reconstruction
	 * @param recImage Matrix with the reconstructed data
	 * @param objPixSizes Pixel sizes in X and Y directions
	 * @param convergenceCurve The convergence curve if available
	 */
	Reconst(std::string scanID,
			Eigen::MatrixXd recImage,
			const std::array<double, 2>& objPixSizes,
			std::vector<double> convergenceCurve=std::vector<double>(0)): Object2D(recImage, objPixSizes),
										   scanID(scanID),
										   convergenceCurve(convergenceCurve)
    {
    }

	std::vector<double> getConvergenceCurve() const;

private:
	std::string scanID; /** Identifier of the scan */
	std::vector<double> convergenceCurve; /** Curve of the convergence in case of iterative reconstruction if available */
};



