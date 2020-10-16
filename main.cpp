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

	Gen1CT ct(300, 300);
	ct.putObject(&phantom);

	constexpr int numProjections{100};
	std::vector<double> angles(180);
	for(int i=0; i<numProjections; i++){angles[i]=i*180/M_PI;}
	//ct.measure(angles);

	int tmpi;
	std::cin>>tmpi;
	return 0;
}
