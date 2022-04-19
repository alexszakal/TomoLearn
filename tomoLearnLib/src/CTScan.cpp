#include <CTScan.hpp>
#include <Object2D.hpp>

#include <iostream>

/**
 * Default CTScan constructor.
 */
CTScan::CTScan():Object2D(),scanID("Empty"),angles(),I0(0){
}

/**
 * CTScan constructor to initialize a CTScan object from an Eigen matrix data source holding the sinogram data
 * @param scanID  Identifier ID of the scan
 * @param sinogram Sinogram data
 * @param detWidth Width of the detector array [mm]
 * @param angles Eigen::Vector holding the projection angle [rad] values
 * @param I0 The source strength used for the scan
 */
CTScan::CTScan(std::string scanID,
		       Eigen::MatrixXd sinogram,
			   double detWidth,
			   const Eigen::VectorXd& angles,
			   double I0):
				        	                 Object2D(sinogram, detWidth, angles),
											 scanID{scanID},
											 angles{angles},
											 I0(I0)
{

}

/**
 * CTScan constructor which calculates the sum of analytical sinograms of ellises.
 * @param scanID Identifier ID of the scan
 * @param detWidth Width of the detector array [mm]
 * @param numDetPixels Number of pixels of the detector array
 * @param angles Eigen::Vector holding the projection angle [rad] values
 * @param rhos Attenuation coeff. of the ellipses [1/mm]
 * @param alphas Inclination angle of ellipses [rad]
 * @param centers Centers of ellipses [mm]
 * @param axes Lengths of the ellipses axes [mm]
 * @param I0 The source strength used for the scan
 */
CTScan::CTScan(std::string scanID,
			double detWidth,
			int numDetPixels,
			const Eigen::VectorXd& angles,
			const std::vector<double>& rhos,
			const std::vector<double>& alphas,
			const std::vector< std::array<double,2> >& centers,
			const std::vector< std::array<double,2> >& axes,
			double I0):Object2D(Eigen::MatrixXd::Zero(numDetPixels, angles.size()),
					                                                   detWidth,
																	   angles),
														      scanID{scanID},
															  angles{angles},
															  I0(I0){

	const Eigen::MatrixXd& objDataRef = getDataAsEigenMatrixRef();

	for(int ellipseIdx =0; static_cast<size_t>(ellipseIdx) < rhos.size(); ++ellipseIdx){
		double x0=centers[static_cast<size_t>(ellipseIdx)][0];
		double y0=centers[static_cast<size_t>(ellipseIdx)][1];
		double A=axes[static_cast<size_t>(ellipseIdx)][0];
		double B=axes[static_cast<size_t>(ellipseIdx)][1];

		double s = std::sqrt( std::pow(x0,2) + std::pow(y0,2) );
		double gamma{0.0};
		if( (x0 != 0.0) || (y0 != 0.0) )
			gamma = std:: atan2(y0, x0);

		for(int angI=0; angI<angles.size(); ++angI){
			double thetaMinusAlpha = angles[static_cast<long>(angI)] - alphas[static_cast<size_t>(ellipseIdx)];
			double aSquaredTheta = std::pow(A*std::cos(thetaMinusAlpha), 2) + std::pow(B*std::sin(thetaMinusAlpha), 2) ;
			for(int pixI=0; pixI<numDetPixels; ++pixI){
				double t=getXValueAtPix(pixI) - s* std::cos(gamma - getYValueAtPix(angI));

				if(t*t <= aSquaredTheta){
					setData(pixI, angI, objDataRef(pixI, angI) + 2*rhos[static_cast<size_t>(ellipseIdx)]*A*B/aSquaredTheta*std::sqrt(aSquaredTheta - std::pow(t,2)) );
				}

			}
		}
	}
}

/**
 * Returns a const reference to the vector holding the projection angle values
 * @return Const reference holding the projection angle values
 */
const Eigen::VectorXd& CTScan::getAnglesConstRef() const{
	return angles;
}

/**
 * Getter function of the detector width
 * @return Detector width in [mm]
 */
double CTScan::getDetWidth() const{
	return getNumberOfPixels()[0] * getPixSizes()[0];
}

/**
 * Getter function of the source strength used for the scan
 * @return Source strength [1]
 */
double CTScan::getI0() const{
	return I0;
}

/**
 * Elementwise add two CTScan objects
 * @param lhs CTScan object on the left handside
 * @param rhs CTScan object on the right handside
 * @return CTSCAN object holding the sum
 */
CTScan operator+(const CTScan& lhs, double rhs){
	return CTScan( lhs.scanID,
		       lhs.getDataAsEigenMatrixRef().array() + rhs,
			   lhs.getDetWidth(),
			   lhs.getAnglesConstRef(),
			   lhs.I0);
}

/**
 * Elementwise division of two CTScan objects
 * @param lhs CTScan object on the left handside
 * @param rhs CTScan object on the right handside
 * @return CTSCAN object holding the result of division
 */
CTScan operator/(const CTScan& lhs, const CTScan& rhs){
	return CTScan( lhs.scanID,
		       lhs.getDataAsEigenMatrixRef().array() / (rhs.getDataAsEigenMatrixRef().array()),
			   lhs.getDetWidth(),
			   lhs.getAnglesConstRef(),
			   lhs.I0);
}

/**
 * Elementwise subtract two CTScan objects
 * @param lhs CTScan object on the left handside
 * @param rhs CTScan object on the right handside
 * @return CTSCAN object holding the difference object
 */
CTScan operator-(const CTScan& lhs, const CTScan& rhs){
	return CTScan( lhs.scanID,
		       lhs.getDataAsEigenMatrixRef().array() - (rhs.getDataAsEigenMatrixRef().array()),
			   lhs.getDetWidth(),
			   lhs.getAnglesConstRef(),
			   lhs.I0);
}

/**
 * Convert the data held in CTScan from counts to line integral values.
 *
 * Converts the values IN PLACE!!!!!
 * Calculates the (-1)*std::log(CTScan(i,j)/I0) value
 */
void CTScan::convertToLineIntegrals(){

	const Eigen::MatrixXd& objDataRef = getDataAsEigenMatrixRef();

	//Calculate the line integrals from the measured counts
	for(int i=0; i<objDataRef.rows(); i++){
		for(int j=0; j<objDataRef.cols(); j++){
			double lineIntegralValue = (-1)* std::log( objDataRef(i,j) / I0 );
			if (lineIntegralValue<0.0){    //The line integral has to be non-negative. It can become neg. because of statistics
				setData(i,j, 0);
			} else{
				setData(i,j, lineIntegralValue);
			}
		}
	}
}

/**
 * Elementwise multiply two CTScans
 * @param lhs Left handside operand of multiply
 * @param rhs Right handside operand of multiply
 * @return Elementwise multiplied CTScan object
 */
CTScan operator*(const CTScan& lhs, const CTScan& rhs){
	return CTScan( lhs.scanID,
				   lhs.getDataAsEigenMatrixRef().array() * rhs.getDataAsEigenMatrixRef().array(),
			       lhs.getDetWidth(),
			       lhs.getAnglesConstRef(),
			       lhs.I0);
}

/**
 * Multiply a CTScan object with a double
 * @param lhs multiplier
 * @param rhs CTScan object to be multiplied
 * @return Multiplied CTScan object
 */
CTScan operator*(double lhs, const CTScan& rhs){
	return CTScan( rhs.scanID,
			       rhs.getDataAsEigenMatrixRef().array()*lhs,
				   rhs.getDetWidth(),
				   rhs.getAnglesConstRef(),
				   rhs.I0);
}

/**
 * Calculate elementvise exp(CTScan)
 * @return CTScan object with the exponentiated elements
 */
CTScan CTScan::exp(){
	return CTScan(scanID,
			      getDataAsEigenMatrixRef().array().exp(),
				  getDetWidth(),
				  angles,
				  I0);
}
































