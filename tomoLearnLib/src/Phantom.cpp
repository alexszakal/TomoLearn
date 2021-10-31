#include <Phantom.hpp>
#include <iostream>
#include <cmath>

//Phantom::Phantom():Object2D(){
//
//}

Phantom::Phantom():Object2D(),label("Empty"){
}

Phantom::Phantom(const std::string& label,
		         const std::string& imageFilePath,
				 const std::array<double, 2>& objPixSizes):
				                                          Object2D{imageFilePath, objPixSizes},
														  label{label}{
}

Phantom::Phantom(const std::string& label,
		         const Eigen::MatrixXd inData,
				 const std::array<double, 2> objPixSizes):
						                                  Object2D{inData, objPixSizes},
														  label{label}{
}

Phantom::Phantom(std::string label, const Object2D& dataPar):Object2D{dataPar},
		                                                     label{label}{
}

Phantom::Phantom(std::string label,
		const std::array<int,2>& numberOfPixels,
		const std::array<double,2>& objPixSizes,
		const std::vector<double>& rhos,
		const std::vector<double>& alphas,
		const std::vector< std::array<double,2> >& centers,
        const std::vector< std::array<double,2> >& axes ):Object2D{numberOfPixels, objPixSizes},
        		                                          label{label}{

        const Eigen::MatrixXd& objDataRef = getDataAsEigenMatrixRef();

        for(int ellipseIdx =0; ellipseIdx<rhos.size(); ++ellipseIdx){
        	for(int i=0; i<numberOfPixels[0]; ++i){
        		for(int j=0; j<numberOfPixels[1]; ++j){
        			double xEllipse=   (getXValueAtPix(i)-centers[ellipseIdx][0])*cos(alphas[ellipseIdx]) + (getYValueAtPix(j)-centers[ellipseIdx][1])*sin(alphas[ellipseIdx]);
        			double yEllipse=-1*(getXValueAtPix(i)-centers[ellipseIdx][0])*sin(alphas[ellipseIdx]) + (getYValueAtPix(j)-centers[ellipseIdx][1])*cos(alphas[ellipseIdx]);

        			double eVal = pow(xEllipse/axes[ellipseIdx][0], 2) +
        				      pow(yEllipse/axes[ellipseIdx][1], 2);

        			if( eVal <1 ){
        				setData(i,j, objDataRef(i,j)+rhos[ellipseIdx]);
        			}
        		}
        	}
        }



}

std::string Phantom::getLabel() const{
	return label;
}

Phantom Phantom::operator*(double coeff) const {
	return Phantom{label, static_cast<Object2D>(*this) * coeff};
}

Phantom Phantom::operator+(double addVal) const {
	return Phantom{label, static_cast<Object2D>(*this) + addVal};
}

Phantom Phantom::operator-(double subVal) const {
	return *this + (-1.0*subVal);
}








