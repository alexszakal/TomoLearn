#include <iostream>
#include <cmath>
#include <chrono>

#include <CImg.h>
#ifdef Success       //Because otherwise Eigen not compile (EIgen <-> CImg interference)
  #undef Success
#endif

#include <TomoLearnS/Object2D.hpp>
#include <TomoLearnS/Gen1CT.hpp>

//Parallel geometry
int main(){
	std::cout << "Parallel beam projection simulation" << std::endl;
	//Reading Shepp-Logan phantom
	Object2D phantom(std::string("Phantoms/SheppLogan.png") );
	phantom.display("Shepp-Logan phantom");

	Gen1CT ct(95, 128);  //width[mm], pixNum
	ct.putObject(&phantom);

	const int numProjections{110};
	std::vector<double> angles(numProjections);
	for(int i=0; i<numProjections; i++){angles[i]=i/static_cast<double>(numProjections)*M_PI;}
	ct.measure(angles);
	ct.displayMeasurement();

	ct.FBP(std::vector<int>{1024,1024}, std::vector<double>{0.1, 0.1});
	//ct.displayMeasurement();

	ct.backProject(std::vector<int>{1024,1024}, std::vector<double>{0.1, 0.1});

	int tmpi;
	std::cin>>tmpi;
	return 0;
}
