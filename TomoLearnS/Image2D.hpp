#include <CImg.h>

class Image2D{
public:
	Image2D(); //Default constructor, initialize an empty image
	Image2D(std::string imageFilePath);
	void display(std::string title);
private:
	cimg_library::CImg<uint16_t> cimg_image;
	cimg_library::CImgDisplay cimg_window;
};
