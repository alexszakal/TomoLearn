#pragma once
#include <Object2D.hpp>
#include <string>
#include <thread>
#include <array>
#include <Eigen/Dense>

#include <iostream>

/**
 * @brief CTScan class stores the integrated attenuation coefficients or simulated counts
 *
 * If CTScan.I0 == 0 it holds the numerically integrated values of the line integrals
 * If CTScan.I0 != 0 the value of the line integrals is computed by applying Poisson
 *                   statistics to model the effects of noise on the reconstructed image.
 */
class CTScan : public Object2D{
public:
	CTScan();
	CTScan(std::string scanID, Eigen::MatrixXd sinogram, double detWidth, const Eigen::VectorXd& angles, double I0);
	CTScan(std::string scanID,
			double detWidth,
			int numDetPixels,
			const Eigen::VectorXd& angles,
			const std::vector<double>& rhos,
			const std::vector<double>& alphas,
			const std::vector< std::array<double,2> >& centers,
			const std::vector< std::array<double,2> >& axes,
			double I0);

	const Eigen::VectorXd& getAnglesConstRef() const;
	double getDetWidth() const;
	double getI0() const;
	void convertToLineIntegrals();

	friend CTScan operator/(const CTScan& lhs, const CTScan& rhs);
	friend CTScan operator-(const CTScan& lhs, const CTScan& rhs);
	friend CTScan operator+(const CTScan& lhs, double rhs);
	friend CTScan operator*(const CTScan& lhs, const CTScan& rhs);
	friend CTScan operator*(double lhs, const CTScan& rhs);
	CTScan exp();

private:
	std::string scanID;   /** Unique ID of the CTScan */
	Eigen::VectorXd angles;  /** Projection angles of the CTScan [rad] */
	double I0;  /** Source strength used for the CTScan */
};
