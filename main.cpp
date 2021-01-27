#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>

#include <CImg.h>
#ifdef Success       //Because otherwise Eigen not compile (EIgen <-> CImg interference)
  #undef Success
#endif

#include <TomoLearnS/Object2D.hpp>
#include <TomoLearnS/Gen1CT.hpp>

#include <matplotlibcpp/matplotlibcpp_old.h>

void testRadonTransform();
void testFBP();

//Parallel geometry
int main(){
	//testRadonTransform();

    testFBP();

	int tmpi;
	std::cin>>tmpi;

	return 0;
}

void testFBP(){
	/**
	 * Test the Filtered Backprojection algorithm with a Shepp-Logan phantom
	 */

	std::cout << "Parallel beam FBP simulation" << std::endl;
	//Reading Shepp-Logan phantom
	//Object2D phantom(std::string("Phantoms/ModifiedSheppLogan_asymmetric.png") );
	//phantom.display("Modified Shepp-Logan phantom");

	Gen1CT ct(110, 1024);  //width[mm], pixNum
	//ct.addPhantom("SL", "Phantoms/SheppLogan.png");
	//ct.addPhantom("SL", "Phantoms/SheppLogan_asymmetric.png");
	//ct.addPhantom("SL", "Phantoms/ModifiedSheppLogan.png");
	ct.addPhantom("SL", "Phantoms/ModifiedSheppLogan_asymmetric.png");

	ct.displayPhantom("SL");

	const int numProjections{180};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0, 180/180*M_PI);

	//ct.measure("SL", angles, "SLsinogram");
/*	ct.displayMeasurement("SLsinogram");

	ct.FBP(std::vector<int>{1024,1024}, std::vector<double>{0.1, 0.1});
	//ct.displayMeasurement();

	ct.backProject(std::vector<int>{1024,1024}, std::vector<double>{0.1, 0.1});
*/

	int tmpi;
	std::cin>>tmpi;

}

void testRadonTransform(){
	/**
	 * Compare the numerical and analytic Radon transform of an ellipse
	 */
/*
	std::cout << "Parallel beam projection simulation" << std::endl;
	//Reading single ellipse phantom
	Object2D phantom(std::string("Phantoms/singleEllipse.png") );
	phantom.display("Single Ellipse phantom");

	//generate RadonTransform
	Gen1CT ct(95, 128);  //width[mm], pixNum
	ct.addPhantom(&phantom);
	const int numProjections{180};
	std::vector<double> angles(numProjections);
	for(int i=0; i<numProjections; i++){angles[i]=i/static_cast<double>(numProjections)*M_PI;}
	ct.measure(angles);
	ct.displayMeasurement();

	//Calculate the Radon analytically
	//std::array<int,2> numberOfPixels = phantom.getNumberOfPixels();
	auto numberOfPixels = phantom.getNumberOfPixels();
	auto pixSizes = phantom.getPixSizes();
	double x0=0.2;  // 1.0 = edge of image (MATLAB phantom convention)
	double y0=0.3;  //    --||--
	double a=0.3;	//relative to the full image size (MATLAB phantom convention)
	double b=0.45;	//     --||--
	double A=255.0;
	double phi=30.0 /180*M_PI;
	double theta=5.0 /180*M_PI;

	//Search the index of the closest measured angle to theta
	int angleIndex = std::distance(angles.begin(), std::min_element(angles.begin(), angles.end(), [theta](double a, double b){return std::abs(a-theta) < std::abs(b-theta); } ) );
	std::cout << std::endl << "AngleIndex: " << angleIndex << " angle=" << angles[angleIndex];

	Eigen::VectorXd projection = ct.getSinogram().col(angleIndex);
	Eigen::VectorXd aProjection = projection*0; //analytical projection
	std::vector<double> pixPositions = ct.getPixPositions();

	double aParam = std::sqrt(  std::pow((a/2*numberOfPixels[0]*pixSizes[0]) * std::cos(theta-phi),2) + std::pow( (b/2*numberOfPixels[1]*pixSizes[1]) * std::sin(theta-phi),2) );
	double s = std::sqrt( std::pow(x0/2*numberOfPixels[0]*pixSizes[0], 2) + std::pow(y0/2*numberOfPixels[1]*pixSizes[1], 2) );

	double gamma{0.0};
	if( (x0 != 0.0) || (y0 != 0.0) )
		gamma = std::atan2(y0, x0);

	std::cout << std::endl << "a=" << a*numberOfPixels[0]*pixSizes[0]/2;
	std::cout << std::endl << "b=" << b*numberOfPixels[1]*pixSizes[1]/2;

	for(int i=0; i<projection.rows(); ++i){
		if (std::fabs( pixPositions[i] - s*cos(gamma-theta) ) < aParam){
			aProjection(i)=2*A*(a*numberOfPixels[0]*pixSizes[0]/2) * (b*numberOfPixels[1]*pixSizes[1]/2) / std::pow(aParam,2) *std::sqrt(std::pow(aParam,2) - std::pow(pixPositions[i] - s* std::cos(gamma-theta),2 ));
		}
		std::cout << std::endl << "numeric: " << projection(i) << "  analytic: " << aProjection(i);
	}

	double coeff =static_cast<double>(projection.maxCoeff())/aProjection.maxCoeff();
	std::cout << std::endl << "Coefficient max(numerical)/max(analytical): " << coeff;

	matplotlibcpp::figure(1);
	matplotlibcpp::plot(std::vector<float> (&projection[0],  projection.data()+  projection.cols()*projection.rows()) );
	matplotlibcpp::plot(std::vector<float> (&aProjection[0], aProjection.data()+aProjection.cols()*aProjection.rows()) );
	matplotlibcpp::show(False);

	int tmpi;
	std::cin>>tmpi;

*/
}


