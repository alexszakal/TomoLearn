#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>

#include <CImg.h>

#ifdef Success       //Because otherwise Eigen not compile (EIgen <-> CImg interference)
  #undef Success
#endif

#include <matplot/matplot.h>

#include <Object2D.hpp>
#include <Gen1CT.hpp>
#include <Phantom.hpp>

#include <config.h>

#include <matplot/matplot.h>


void testSPS(const std::string& phantomName, projectorType measureAlgo, projectorType projectAlgo,
		  backprojectorType backprojectAlgo);

int main(){

#if ENABLE_CUDA
	std::cout << "\n \n CUDA enabled!!!!" ;
#else
	std::cout << "\n \n CUDA disabled!!!" ;
#endif

	//rayDriven Projector and pixelDriven BackProjector Standard config.
	//testSPS("modSL_symm", projectorType::rayDriven, projectorType::rayDriven, backprojectorType::pixelDriven);

	#if ENABLE_CUDA
		//rayDriven Projector on GPU and pixelDriven BackProjector Standard config.
		//testSPS("modSL_symm", projectorType::rayDriven, projectorType::rayDriven_GPU, backprojectorType::pixelDriven);

		//rayDriven Projector and Backprojector on GPU
		testSPS("modSL_symm", projectorType::rayDriven, projectorType::rayDriven_GPU, backprojectorType::rayDriven_GPU);
	#endif

	//testSPS("modSL_symm", projectorType::rayDriven, projectorType::rayDriven, backprojectorType::rayDriven);

	//testSPS("modSL_symm", projectorType::rayDriven, projectorType::pixelDriven, backprojectorType::pixelDriven);

	//testSPS("modSL_symm", projectorType::rayDriven, projectorType::pixelDriven, backprojectorType::rayDriven);

	//testSPS("modSL_symm", projectorType::rayDriven, projectorType::Siddon, backprojectorType::pixelDriven);

	//testSPS("modSL_symm", projectorType::rayDriven, projectorType::Siddon, backprojectorType::rayDriven);

	std::cout<<"\n Press ENTER to finish";
	std::cin.get();

	return 0;
}

void testSPS(const std::string& phantomName,
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
	//ct.setI0(0.0);

	ct.measure(phantomName, angles, "Sinogram", measureAlgo);

	ct.displayMeasurement("Sinogram");

	std::cout<<"\nStart SPS reconstruction";
	ct.SPSReconst("Sinogram", std::array<int, 2> { 512, 512}, // 1024 x 1024 pixel, 0.1mm felbontas
			std::array<double, 2> { 0.199, 0.199}, projectAlgo, backprojectAlgo, "RecImage", 110,
			regularizerType::Huber, 1000, 0.004, phantomName);

	ct.Gen1CT::displayReconstruction("RecImage");

	ct.compareRowPhantomAndReconst('Y', -31.0, phantomName, "RecImage");

	auto h=matplot::figure(true);
	auto p1 = matplot::plot( ct.getConvergenceCurve("RecImage"), "-og" );
	matplot::title("Convergence of SPSreconstruction");
	matplot::xlabel("Iteration number [1]");
	matplot::ylabel("Difference L2 norm [a.u.]");
	h->show();


}
