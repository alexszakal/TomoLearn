#include <TomoLearnS/Object2D.hpp>

Object2D::Object2D(std::string imageFilePath, std::array<double, 2> objPixSizes):objPixSizes{objPixSizes}{
	cimg_image=cimg_library::CImg<uint16_t>(imageFilePath.c_str());
}

void Object2D::display(std::string title){
	cimg_window=cimg_library::CImgDisplay(cimg_image, title.c_str());
	cimg_window.wait();
}
