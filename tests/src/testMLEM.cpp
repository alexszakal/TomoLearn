#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>

#include <CImg.h>

#ifdef Success       //Because otherwise Eigen not compile (EIgen <-> CImg interference)
  #undef Success
#endif

#include <Object2D.hpp>
#include <Gen1CT.hpp>

#include <matplotlibcpp.h>

#include <config.h>


void testMLEM(const std::string& phantomName, projectorType measureAlgo, projectorType projectAlgo,
		  const std::string& backprojectAlgo);

int main(){

#if ENABLE_CUDA
	std::cout << "\n \n CUDA enabled!!!!" ;
#else
	std::cout << "\n \n CUDA disabled!!!" ;
#endif

	//Test the naive implementations with interpolations
//	testFBP("modSL_symm", "withInterpolation", "backProject_interpol");

//	testFBP("modSL_symm", "haoGaoProject", "backProject_interpol");

//	testFBP("modSL_symm", "withInterpolation", "backProject_HaoGao_CPU");

//	testMLEM("modSL_symm", "haoGaoProject");
	testMLEM("modSL_symm", projectorType::pixelDriven, projectorType::pixelDriven, "backprojectSimple");

	std::cin.ignore();

	return 0;
}

void testMLEM(const std::string& phantomName,
		      projectorType measureAlgo,
			  projectorType projectAlgo,
			  const std::string& backprojectAlgo){
	/**
	 * Test the MLEM algorithm with a given phantom
	 */

	std::cout << "MLEM reconstruction simulation" << std::endl;

	int detWidthInMM { 150 };
	int detPixNum { 1024 };
	Gen1CT ct(detWidthInMM, detPixNum);

	//Reading Shepp-Logan phantom
	ct.addPhantom("SL", "Phantoms/SheppLogan_HU.png");
	ct.addPhantom("SL_asym", "Phantoms/SheppLogan_asymmetric_HU.png");
	ct.addPhantom("modSL_symm", "Phantoms/ModifiedSheppLogan_HU.png");
	ct.addPhantom("modSL_asym", "Phantoms/ModifiedSheppLogan_asymmetric_HU.png"); //default pixSize: 0.1mm x 0.1mm
	ct.addPhantom("SD", "Phantoms/SingleDot.png"); //Single dot Phantom
	ct.addPhantom("rectangle", "Phantoms/rectangle.png", {0.025,0.025}); //Single rectangle with 400HU CT number in  the center

	ct.displayPhantom(phantomName);

	const int numProjections{180};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0.0/180.0 * M_PI,
			(1.0 - 1.0/numProjections) * M_PI);

	ct.setI0(8e4);
	ct.setI0(0.0);

	ct.measure(phantomName, angles, "Sinogram", measureAlgo);

	ct.displayMeasurement("Sinogram");

	ct.MLEMReconst("Sinogram", std::array<int, 2> { 256, 256}, //jo 256 x 256 pixel, 0.4 felbontas
			std::array<double, 2> { 0.4, 0.4 }, projectAlgo, backprojectAlgo, "RecImage", 100);

	ct.Gen1CT::displayReconstruction("RecImage");

	ct.compareRowPhantomAndReconst('Y', -31.0, phantomName, "RecImage");

	std::cout<<"\n Press ENTER to continue";
	std::cin.get();
}
