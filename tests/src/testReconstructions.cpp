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

int main(){

	#if ENABLE_CUDA
		std::cout << "\n \n CUDA enabled!!!!" ;
	#else
		std::cout << "\n \n CUDA disabled!!!" ;
	#endif

	//Test configuration
	std::string phantomFileName = "Phantoms/ModifiedSheppLogan_HU.png";
	projectorType measureAlgo = projectorType::pixelDriven;
	projectorType projectAlgo = projectorType::rayDriven;
	backprojectorType backprojectAlgo = backprojectorType::rayDriven;
	double I0 = 3e2;

	int MLEM_numIterations=100;
	int SPS_numIterations=100;
	//TEST CONFIG END

	//Generate Sinogram
	double detWidthInMM { 145.3 };
	int detPixNum { 1453 };
	Gen1CT ct(detWidthInMM, detPixNum);

	ct.addPhantom("Phantom", phantomFileName, {0.1, 0.1}, true);
	ct.displayPhantom("Phantom");

	const int numProjections{180*2};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0.0/180.0 * M_PI,
				(1.0 - 1.0/numProjections) * M_PI);

	ct.setI0(I0);

	ct.measure("Phantom", angles, "Sinogram", measureAlgo);

	ct.displayMeasurement("Sinogram");

	//Reconstruct with FBP
	ct.filteredBackProject("Sinogram", std::array<int, 2> { 512, 512}, //jo 256 x 256 pixel, 0.4 felbontas
				std::array<double, 2> { 0.2, 0.2 }, FilterType::Hann, 0.5, backprojectAlgo,
				"FBP_RecImage", "Phantom");
	ct.displayReconstruction("FBP_RecImage");
	auto h=matplot::figure(true);
	auto p0 = matplot::plot( std::vector<double>(std::max(MLEM_numIterations, SPS_numIterations), ct.getConvergenceCurve("FBP_RecImage")[0]), "-xr" );
	matplot::hold(true);
	matplot::title("Convergence of different reconstructions");
	matplot::xlabel("Iteration number [1]");
	matplot::ylabel("Difference L2 norm [a.u.]");
	p0->display_name("FBP value");
	matplot::legend();
	h->draw();

	//Reconstruct with MLEM
	ct.MLEMReconst("Sinogram", std::array<int, 2> { 512, 512}, // 1024 x 1024 pixel, 0.1mm felbontas
				std::array<double, 2> { 0.2, 0.2}, projectAlgo, backprojectAlgo, "MLEM_RecImage", MLEM_numIterations, "Phantom");
	ct.displayReconstruction("MLEM_RecImage");


	auto p1 = matplot::plot( ct.getConvergenceCurve("MLEM_RecImage"), "-og" );
	p1->display_name("MLEM convergence");
	h->draw();

	//Reconstuct with SPS
	ct.SPSReconst("Sinogram", std::array<int, 2> { 512, 512}, // 1024 x 1024 pixel, 0.1mm felbontas
						std::array<double, 2> { 0.199, 0.199}, projectAlgo, backprojectAlgo, "SPS_RecImage", SPS_numIterations,
						regularizerType::Huber, 3000, 0.004, "Phantom");
	ct.displayReconstruction("SPS_RecImage");

	auto p2 = matplot::plot( ct.getConvergenceCurve("SPS_RecImage"), "-ob" );
	p2->display_name("SPS convergence");
	h->draw();

	std::cin.get();

}

