
#include <Gen1CT.hpp>

#include <cstring>


int main(int argc, char* argv[]){

	if(argc!=2){
		std::cout << "\nWrong number of parameters!! Aborting.";
		return 1;
	}

	std::cout<<"\nprogram: " << argv[0];

	double detWidthInMM { 145.3 };
	int detPixNum { 1453 };
	Gen1CT ct(detWidthInMM, detPixNum);

	ct.addPhantom("modSL_symm", "Phantoms/ModifiedSheppLogan_HU.png", {0.1, 0.1}, true);

    //ct.displayPhantom("modSL_symm");

	const int numProjections{180*2};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0.0/180.0 * M_PI,
				(1.0 - 1.0/numProjections) * M_PI);

	//ct.setI0(1e5);
	ct.setI0(0.0);

	if(strcmp(argv[1],"1") == 0){
		ct.measure("modSL_symm", angles, "Sinogram", projectorType::rayDriven);
	}
	else if(strcmp(argv[1], "2")==0){
		ct.measure("modSL_symm", angles, "Sinogram", projectorType::rayDrivenOptimized);
	}
	else if(strcmp(argv[1], "3")==0){
		ct.displayPhantom("modSL_symm");

		ct.measure("modSL_symm", angles, "Sinogram", projectorType::rayDriven);
		ct.displayMeasurement("Sinogram");
		CTScan rayDriven = ct.getMeasurement("Sinogram");

		ct.measure("modSL_symm", angles, "SinogramOptimized", projectorType::rayDrivenOptimized);
		ct.displayMeasurement("SinogramOptimized");
		CTScan rayDrivenOptimized = ct.getMeasurement("SinogramOptimized");

		CTScan diffImage = rayDriven - rayDrivenOptimized;
		diffImage.display("diffImage");
	}
	else{
		std::cout<<"\nWrong parameter!!";
	}

	//ct.displayMeasurement("Sinogram");

}
