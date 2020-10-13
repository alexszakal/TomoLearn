#include <iostream>
#include <CImg.h>
#include <TomoLearnS/Image2D.hpp>

//Parallel geometry
int main(){
	std::cout << "Parallel geometry simulation" << std::endl;

	//Reading Shepp-Logan phantom
	Image2D phantom("Phantoms/SheppLogan.png");
	phantom.display("Shepp-Logan phantom");



	return 0;
}
