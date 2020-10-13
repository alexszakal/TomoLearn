#include <iostream>
#include <CImg.h>
#include <TomoLearnS/Object2D.hpp>

//Parallel geometry
int main(){
	std::cout << "Parallel beam projection simulation" << std::endl;
	//Reading Shepp-Logan phantom
	Object2D phantom("Phantoms/SheppLogan.png");
	phantom.display("Shepp-Logan phantom");

	return 0;
}
