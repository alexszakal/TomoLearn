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

void testFBP(const std::string& phantomName, projectorType projectAlgo, backprojectorType backprojectAlgo);

//TODO: A ct.compareRowPhantomAndReconst() Mukodjon. HA fajlbol olvasunk, akkor 1000-et ki kell vonni, mert akkor kapjuk meg HU unitban!

//DEBUG: A szurt szinogram eltunik amikor a visszaallitas megjelenik
//TODO: Valahogy a szurt szinogramokat is el kell menteni (lehetne egy map, ahol a key a filter osztaly?? )

//TODO: Visszavetitest felgyorsitani
//TODO: Gyorsabb elorevetites a cache jobb hasznalataval

//TTOK sanitizers:
//-Dokumentacioba Properties -> C/C++Build -> CMake4eclipse -> Symbols -> Cmake cache entries -> ENABLE_SANITIZER_ADDRESS:BOOL:ON
//- Az address symbolizer is kell ahhoz hogy kodsorokat irjon ki: LAunch config. -> Environmentbe:
                                                           // -> ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer

//Parallel geometry
int main(){
//	testRadonTransform("SL", "swithInterpolation");

#if ENABLE_CUDA
	std::cout << "\n \n CUDA enabled!!!!" ;
#else
	std::cout << "\n \n CUDA disabled!!!" ;
#endif

	//Works: rayDriven Projector and pixelDriven BackProjector
	testFBP("modSL_symm", projectorType::rayDriven, backprojectorType::pixelDriven);

	//Works: rayDriven Projector and rayDriven BackProjector
	//testFBP("modSL_symm", projectorType::rayDriven, backprojectorType::rayDriven);

	//Works: pixelDriven Projector and pixelDriven BackProjector
	//testFBP("modSL_symm", projectorType::pixelDriven, backprojectorType::pixelDriven);

	//Works: pixelDriven Projector and rayDriven BackProjector
	//testFBP("modSL_symm", projectorType::pixelDriven, backprojectorType::rayDriven);

	//Works: Siddon Projector and pixelDriven BackProjector
	//testFBP("modSL_symm", projectorType::Siddon, backprojectorType::pixelDriven);

	//Works: Siddon Projector and rayDriven BackProjector
	//testFBP("modSL_symm", projectorType::Siddon, backprojectorType::rayDriven);

	std::cin.ignore();

	return 0;
}

void testFBP(const std::string& phantomName,
		     projectorType projectAlgo,
			backprojectorType backprojectAlgo){
	/**
	 * Test the Filtered Backprojection algorithm with a Shepp-Logan phantom
	 */

	std::cout << "Parallel beam FBP simulation" << std::endl;

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

	ct.displayPhantom(phantomName);

	const int numProjections{180*2};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0.0/180.0 * M_PI,
			(1.0 - 1.0/numProjections) * M_PI);

	ct.setI0(1e2);
	//ct.setI0(0.0);

	ct.measure(phantomName, angles, "Sinogram", projectAlgo);

	ct.displayMeasurement("Sinogram");

	ct.filteredBackProject("Sinogram", std::array<int, 2> { 512, 512}, //jo 256 x 256 pixel, 0.4 felbontas
			std::array<double, 2> { 0.2, 0.2 }, FilterType::Hann, 0.5, backprojectAlgo,
			"RecImage");

	ct.Gen1CT::displayReconstruction("RecImage");

	ct.compareRowPhantomAndReconst('Y', -31.0, phantomName, "RecImage");

	std::cout<<"\n Press ENTER to continue";
	std::cin.get();
}



