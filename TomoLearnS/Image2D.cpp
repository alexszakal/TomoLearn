#include <CImg.h>
#include <cstdint>
#include <string>
#include <TomoLearnS/Image2D.hpp>

Image2D::Image2D(){
	cimg_image = cimg_library::CImg<uint16_t>();
}

Image2D::Image2D(std::string imageFilePath){
	cimg_image=cimg_library::CImg<uint16_t>(imageFilePath.c_str());
}

void Image2D::display(std::string title){
	cimg_window=cimg_library::CImgDisplay(cimg_image, title.c_str());
	cimg_window.wait();
}
