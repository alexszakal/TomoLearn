#include <iostream>
#include <cmath>
#include <CImg.h>
#include <TomoLearnS/Object2D.hpp>
#include <TomoLearnS/Gen1CT.hpp>

//Parallel geometry
int main(){
	std::cout << "Parallel beam projection simulation" << std::endl;
	//Reading Shepp-Logan phantom
	Object2D phantom(std::string("Phantoms/SheppLogan.png") );
	phantom.display("Shepp-Logan phantom");

	Gen1CT ct(10, 10);
	ct.putObject(&phantom);

	constexpr int numProjections{10};
	std::vector<double> angles(numProjections);
	for(int i=0; i<numProjections; i++){angles[i]=(i+20.0)/180.0*M_PI;}
	ct.measure(angles);

	int tmpi;
	std::cin>>tmpi;
	return 0;
}
