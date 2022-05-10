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
#include <Phantom.hpp>

#include <matplotlibcpp.h>

#include <config.h>


void testMLEM(const std::string& phantomName, projectorType measureAlgo, projectorType projectAlgo,
		  backprojectorType backprojectAlgo);

int main(){

#if ENABLE_CUDA
	std::cout << "\n \n CUDA enabled!!!!" ;
#else
	std::cout << "\n \n CUDA disabled!!!" ;
#endif

	//Works: rayDriven Projector and pixelDriven BackProjector Standard config.
	testMLEM("modSL_symm", projectorType::rayDriven, projectorType::rayDriven, backprojectorType::pixelDriven);

	//Works
	//testMLEM("modSL_symm", projectorType::rayDriven, projectorType::rayDriven, backprojectorType::rayDriven);

	//Works
	//testMLEM("modSL_symm", projectorType::rayDriven, projectorType::pixelDriven, backprojectorType::pixelDriven);

	//Works
	//testMLEM("modSL_symm", projectorType::rayDriven, projectorType::pixelDriven, backprojectorType::rayDriven);

	//Works
	//testMLEM("modSL_symm", projectorType::rayDriven, projectorType::Siddon, backprojectorType::pixelDriven);

	//Works
	//testMLEM("modSL_symm", projectorType::rayDriven, projectorType::Siddon, backprojectorType::rayDriven);

	std::cin.ignore();

	return 0;
}

void testMLEM(const std::string& phantomName,
		      projectorType measureAlgo,
			  projectorType projectAlgo,
			  backprojectorType backprojectAlgo){
	/**
	 * Test the MLEM algorithm with a given phantom
	 */

	std::cout << "MLEM reconstruction simulation" << std::endl;

	double detWidthInMM { 145.3 };
	int detPixNum { 1453 };
	Gen1CT ct(detWidthInMM, detPixNum);

	//Reading Shepp-Logan phantom
	ct.addPhantom("SL", "Phantoms/SheppLogan_HU.png", {0.1, 0.1}, true);
	ct.addPhantom("SL_asym", "Phantoms/SheppLogan_asymmetric_HU.png", {0.1, 0.1}, true);
	ct.addPhantom("modSL_symm", "Phantoms/ModifiedSheppLogan_HU.png", {0.1, 0.1}, true);
	ct.addPhantom("modSL_asym", "Phantoms/ModifiedSheppLogan_asymmetric_HU.png", {0.1, 0.1}, true); //default pixSize: 0.1mm x 0.1mm
	ct.addPhantom("SD", "Phantoms/SingleDot.png"); //Single dot Phantom
	ct.addPhantom("rectangle", "Phantoms/rectangle.png", {0.025,0.025}); //Single rectangle with 400HU CT number in  the center

	Phantom tmpTest("SL_symm_nonzeroOutside", "Phantoms/ModifiedSheppLogan_HU.png", {0.1, 0.1}, true);
	tmpTest = tmpTest + 0.003;
	ct.addPhantom(tmpTest);

	Phantom centerDotPhantom("centerDotPhantom", { 1024, 1024 }, { 0.1, 0.1 },
			std::vector<double> { 1000 }, //rhos
			std::vector<double> { 0.0 }, //alphas
			std::vector<std::array<double, 2>> { { 0, 0 } },  //centers
			std::vector<std::array<double, 2>> { { 10, 10 } } //axes
			);
	ct.addPhantom(centerDotPhantom);

	ct.displayPhantom(phantomName);

	const int numProjections{180*2};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0.0/180.0 * M_PI,
			(1.0 - 1.0/numProjections) * M_PI);

	ct.setI0(5e2);
	ct.setI0(0.0);

	ct.measure(phantomName, angles, "Sinogram", measureAlgo);

	ct.displayMeasurement("Sinogram");

	ct.MLEMReconst("Sinogram", std::array<int, 2> { 512, 512}, // 1024 x 1024 pixel, 0.1mm felbontas
			std::array<double, 2> { 0.2, 0.2}, projectAlgo, backprojectAlgo, "RecImage", 2);  //optimalis iteracio: ~60

	std::cout << "\n L2 norm: " << ct.compareReconToPhantom("RecImage", phantomName) <<'\n';

	ct.Gen1CT::displayReconstruction("RecImage");

	ct.compareRowPhantomAndReconst('Y', -31.0, phantomName, "RecImage");

	std::cout<<"\n Press ENTER to continue";
	std::cin.get();
}
