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

void testFBP(const std::string& phantomName, const std::string& algoName);

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

	testFBP("modSL_symm", "Siddon");

	std::cin.ignore();

	return 0;
}

void testFBP(const std::string& phantomName, const std::string& algoName){
	/**
	 * Test the Filtered Backprojection algorithm with a Shepp-Logan phantom
	 */

	std::cout << "Parallel beam FBP simulation" << std::endl;

	int detWidthInMM { 110 };
	int detPixNum { 512 };
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
			179.0 / 180 * M_PI);

	ct.setI0(8e4);

	if(algoName == "Siddon"){
		ct.measure_Siddon(phantomName, angles, "Sinogram");
	} else if(algoName == "withInterpolation"){
		ct.measure_withInterpolation(phantomName, angles, "Sinogram");
	} else{
		std::cout << "\nalgoName parameter not recognized. Possible values: \"Siddon\" or \"withInterpolation\" ";
		std::cout << "\nAborting testRadonTransform function";
		return;
	}

	ct.displayMeasurement("Sinogram");

	ct.filteredBackProject("Sinogram", std::array<int, 2> { 1024, 1024 },
			std::array<double, 2> { 0.1, 0.1 }, FilterType::Hann, 0.5,
			"RecImage");
	ct.Gen1CT::displayReconstruction("RecImage");

	//ct.compareRowPhantomAndReconst(821, phantomName, "RecImage");

	int tmpi;
	std::cin>>tmpi;
}



