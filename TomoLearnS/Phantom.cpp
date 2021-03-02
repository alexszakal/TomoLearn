#include <TomoLearnS/Phantom.hpp>
#include <iostream>

//Phantom::Phantom():Object2D(){
//
//}

Phantom::Phantom(const std::string& label,
		         const std::string& imageFilePath,
				 const std::array<double, 2>& objPixSizes):
				                                          Object2D{imageFilePath, objPixSizes},
														  label{label}{

}






