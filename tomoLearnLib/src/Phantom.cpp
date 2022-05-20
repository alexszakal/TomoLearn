#include <Phantom.hpp>
#include <iostream>
#include <cmath>

#include <Eigen/Dense>

Phantom::Phantom():Object2D(),label("Empty"){
}

/***
 * Constructs the Phantom using data from an image file on disk
 * @param label Identifier label of the phantom
 * @param imageFilePath Path to the file to the initializer image
 * @param objPixSizes Size of the pixels in X and Y directions
 * @param convertFromHUtoLA Switch to convert from HU to Linear Attenuation units
 */
Phantom::Phantom(const std::string& label,
		         const std::string& imageFilePath,
				 const std::array<double, 2>& objPixSizes,
				 bool convertFromHUtoLA):
				                                          Object2D{imageFilePath, objPixSizes, convertFromHUtoLA},
														  label{label}{
}

/***
* Constructs the Phantom using data from a matrix
* @param label Identifier label of the phantom
* @param inData Data which is used for initialization
* @param objPixSizes Size of the pixels in X and Y directions
*/
Phantom::Phantom(const std::string& label,
		         const Eigen::MatrixXd inData,
				 const std::array<double, 2> objPixSizes):
						                                  Object2D{inData, objPixSizes},
														  label{label}{
}

/***
 * Constructs the Phantom using data form an Object2D
 * @param label Identifier label of the phantom
 * @param dataPar Object2D type which is used for initialization
 */
Phantom::Phantom(std::string label, const Object2D& dataPar):Object2D{dataPar},
		                                                     label{label}{
}

/***
 * Constructs the phantom as the sum of ellipses
 * @param label Identifier label of the Phantom
 * @param numberOfPixels Number of pixels in X and Y directions
 * @param objPixSizes Pixel size in X and Y direction
 * @param rhos Values inside the ellipses
 * @param alphas Inclination of the ellipses
 * @param centers Positions of ellipe centers
 * @param axes Size of the axes of the ellipses
 */
Phantom::Phantom(std::string label,
		const std::array<int,2>& numberOfPixels,
		const std::array<double,2>& objPixSizes,
		const std::vector<double>& rhos,
		const std::vector<double>& alphas,
		const std::vector< std::array<double,2> >& centers,
        const std::vector< std::array<double,2> >& axes ):Object2D{numberOfPixels, objPixSizes},
        		                                          label{label}{

        const Eigen::MatrixXd& objDataRef = getDataAsEigenMatrixRef();

        for(size_t ellipseIdx =0; ellipseIdx<rhos.size(); ++ellipseIdx){
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

/***
* Calculate the two terms that is required for the quadratic regularization of the SPS algorithm
* @return Array of the two terms, [0] is the term in the numerator [1] is needed in the denominator
*/
std::array<Phantom, 2> Phantom::calculateQuadRegTerms() const{
    auto image = getDataAsEigenMatrixRef();
    auto imageSize = getNumberOfPixels();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> numeratorTerm = Eigen::MatrixXd::Zero(imageSize[0], imageSize[1]);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> denomTerm = Eigen::MatrixXd::Zero(imageSize[0], imageSize[1]);

    //Go through the pixels
    for(int xIdx=0; xIdx<imageSize[0]; ++xIdx){
        for(int yIdx=0; yIdx<imageSize[1]; ++yIdx){
        	for(int xRelIdx = -1; xRelIdx<=1; ++xRelIdx){
        		for(int yRelIdx = -1; yRelIdx<=1; ++yRelIdx){
        		    //calculate the numeratorTerm
        		    if (xRelIdx ==0 and yRelIdx ==0)
        		        continue;
        		    int xNeighborIdx = xIdx+xRelIdx;
        		    int yNeighborIdx = yIdx+yRelIdx;
        		    if (xNeighborIdx < 0 or xNeighborIdx >= imageSize[0] or  //Check if the neighbor is out of the image
        		           yNeighborIdx < 0 or yNeighborIdx >= imageSize[1]){
        		         continue;
        		    }
        		    double omega_k = 1/std::sqrt(xRelIdx*xRelIdx + yRelIdx*yRelIdx);
        		    numeratorTerm(xIdx, yIdx) += omega_k*(image(xIdx, yIdx) - image(xNeighborIdx, yNeighborIdx));

        		    //Calculate the denominatorTerm
        		    denomTerm(xIdx, yIdx) += 2*omega_k;
        		}
        	}
        }
    }

    return std::array<Phantom, 2>{Phantom("numeratorTerm", numeratorTerm),
        		                  Phantom("denomTerm", denomTerm) };
}

/***
 * Calculate the two terms that is required for the regularization with Huber prior in the SPS algorithm
 * @param delta Regularization parameter of Huber function
 * @return Array of the two terms, [0] is the term in the numerator [1] is needed in the denominator
 */
std::array<Phantom, 2> Phantom::calculateHuberRegTerms(double delta) const{
    auto image = getDataAsEigenMatrixRef();
    auto imageSize = getNumberOfPixels();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> numeratorTerm = Eigen::MatrixXd::Zero(imageSize[0], imageSize[1]);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> denomTerm = Eigen::MatrixXd::Zero(imageSize[0], imageSize[1]);

    //Go through the pixels
    for(int xIdx=0; xIdx<imageSize[0]; ++xIdx){
        for(int yIdx=0; yIdx<imageSize[1]; ++yIdx){
        	for(int xRelIdx = -1; xRelIdx<=1; ++xRelIdx){
        		for(int yRelIdx = -1; yRelIdx<=1; ++yRelIdx){
        		    //calculate the numeratorTerm
        		    if (xRelIdx ==0 and yRelIdx ==0)
        		        continue;
        		    int xNeighborIdx = xIdx+xRelIdx;
        		    int yNeighborIdx = yIdx+yRelIdx;
        		    if (xNeighborIdx < 0 or xNeighborIdx >= imageSize[0] or  //Check if the neighbor is out of the image
        		           yNeighborIdx < 0 or yNeighborIdx >= imageSize[1]){
        		         continue;
        		    }
        		    double omega_k = 1/std::sqrt(xRelIdx*xRelIdx + yRelIdx*yRelIdx);
        		    double t = image(xIdx, yIdx) - image(xNeighborIdx, yNeighborIdx);
        		    if ( std::abs(t) <= delta){
        		    	numeratorTerm(xIdx, yIdx) += omega_k * t;
        		    }
        		    else{
        		    	if(t>=0){
        		    		numeratorTerm(xIdx, yIdx) += omega_k * delta;
        		    	}
        		    	else{
        		    		numeratorTerm(xIdx, yIdx) += omega_k * (-1.0) * delta;
        		    	}
        		    }

        		    //Calculate the denominatorTerm
        		    if(std::abs(t) <= delta){
        		    	denomTerm(xIdx, yIdx) += 2*omega_k;
        		    }
        		    else{
        		    	denomTerm(xIdx, yIdx)=0;
        		    }
        		}
        	}
        }
    }

    return std::array<Phantom, 2>{Phantom("numeratorTerm", numeratorTerm),
        		                  Phantom("denomTerm", denomTerm) };
}

std::array<Phantom,2> Phantom::calculateGibbsRegTerms(double delta) const{
    auto image = getDataAsEigenMatrixRef();
    auto imageSize = getNumberOfPixels();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> numeratorTerm = Eigen::MatrixXd::Zero(imageSize[0], imageSize[1]);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> denomTerm = Eigen::MatrixXd::Zero(imageSize[0], imageSize[1]);

    //Go through the pixels
    for(int xIdx=0; xIdx<imageSize[0]; ++xIdx){
        for(int yIdx=0; yIdx<imageSize[1]; ++yIdx){
        	for(int xRelIdx = -1; xRelIdx<=1; ++xRelIdx){
        		for(int yRelIdx = -1; yRelIdx<=1; ++yRelIdx){
        		    //calculate the numeratorTerm
        		    if (xRelIdx ==0 and yRelIdx ==0)
        		        continue;
        		    int xNeighborIdx = xIdx+xRelIdx;
        		    int yNeighborIdx = yIdx+yRelIdx;
        		    if (xNeighborIdx < 0 or xNeighborIdx >= imageSize[0] or  //Check if the neighbor is out of the image
        		           yNeighborIdx < 0 or yNeighborIdx >= imageSize[1]){
        		         continue;
        		    }
        		    double omega_k = 1/std::sqrt(xRelIdx*xRelIdx + yRelIdx*yRelIdx);
        		    double t = image(xIdx, yIdx) - image(xNeighborIdx, yNeighborIdx);

        		    numeratorTerm(xIdx, yIdx) += omega_k * std::abs(delta)*t/(std::abs(delta) + std::abs(t));

        		    //Calculate the denominatorTerm
        		    denomTerm(xIdx, yIdx) += 2*omega_k * std::pow(delta,2) / std::pow(std::abs(delta) + std::abs(t),2);
        		}
        	}
        }
    }

    return std::array<Phantom, 2>{Phantom("numeratorTerm", numeratorTerm),
        		                  Phantom("denomTerm", denomTerm) };
}

/***
 * Get the identifier label of the Phantom
 * @return Identifier label
 */
std::string Phantom::getLabel() const{
	return label;
}

/***
 * Multiply the Phantom with a scalar
 * @param coeff Scalar used for the multiplication
 * @return Copy of the Phantom with the new values
 */
Phantom Phantom::operator*(double coeff) const {
	return Phantom{label, static_cast<Object2D>(*this) * coeff};
}

/***
 * Add a scalar value to the Phantom
 * @param addVal Scalar value to be added
 * @return Copy of the phantom with the new values
 */
Phantom Phantom::operator+(double addVal) const {
	return Phantom{label, static_cast<Object2D>(*this) + addVal};
}

/***
 * Subtract a scalar value from Phantom
 * @param subVal Scalar value to be subtracted
 * @return Copy of the phantom with the new values
 */
Phantom Phantom::operator-(double subVal) const {
	return *this + (-1.0*subVal);
}

/***
 * Divide the elements of two Phantom objects
 * @param lhs Operand on the left handside
 * @param rhs Operand on the right handside
 * @return A Phantom object with the new values
 */
Phantom operator/(const Phantom& lhs, const Phantom& rhs){
	return Phantom(lhs.label,
		         lhs.getDataAsEigenMatrixRef().array() / (rhs.getDataAsEigenMatrixRef().array()),
				 lhs.getPixSizes());
}

/***
 * Multiply the elements of two Phantom objects
 * @param lhs Operand on the left handside
 * @param rhs Operand on the right handside
 * @return A Phantom object with the new values
 */
Phantom operator*(const Phantom& lhs, const Phantom& rhs){
	return Phantom(lhs.label,
		         lhs.getDataAsEigenMatrixRef().array() * (rhs.getDataAsEigenMatrixRef().array()),
				 lhs.getPixSizes());
}

/***
 * Add the elements of two Phantom objects
 * @param lhs Operand on the left handside
 * @param rhs Operand on the right handside
 * @return A Phantom object with the new values
 */
Phantom operator+(const Phantom& lhs, const Phantom& rhs){
	return Phantom(lhs.label,
		         lhs.getDataAsEigenMatrixRef().array() + (rhs.getDataAsEigenMatrixRef().array()),
				 lhs.getPixSizes());
}








