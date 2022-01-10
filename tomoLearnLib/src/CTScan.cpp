#include <CTScan.hpp>
#include <Object2D.hpp>

#include <iostream>

CTScan::CTScan():Object2D(),scanID("Empty"),angles(),I0(0){
}

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
			gamma = std::atan2(y0, x0);

		for(int angI=0; angI<angles.size(); ++angI){
			double aSquaredTheta = std::pow(A*std::cos(angles[angI]), 2) + std::pow(B*std::sin(angles[angI]), 2) ;
			for(int pixI=0; pixI<numDetPixels; ++pixI){
				double t=getXValueAtPix(pixI) - s* std::cos(gamma - getYValueAtPix(angI));

				if(t*t <= aSquaredTheta){
					setData(pixI, angI, objDataRef(pixI, angI) + 2*rhos[ellipseIdx]*A*B/aSquaredTheta*std::sqrt(aSquaredTheta - std::pow(t,2)) );
				}

			}
		}
	}
}

const Eigen::VectorXd& CTScan::getAnglesConstRef() const{
	return angles;
}

double CTScan::getDetWidth() const{
	return getNumberOfPixels()[0] * getPixSizes()[0];
}

double CTScan::getI0() const{
	return I0;
}

CTScan operator+(const CTScan& lhs, double rhs){
	return CTScan( lhs.scanID,
		       lhs.getDataAsEigenMatrixRef().array() + rhs,
			   lhs.getDetWidth(),
			   lhs.getAnglesConstRef(),
			   lhs.I0);
}

CTScan operator/(const CTScan& lhs, const CTScan& rhs){
	return CTScan( lhs.scanID,
		       lhs.getDataAsEigenMatrixRef().array() / (rhs.getDataAsEigenMatrixRef().array()),
			   lhs.getDetWidth(),
			   lhs.getAnglesConstRef(),
			   lhs.I0);
}
