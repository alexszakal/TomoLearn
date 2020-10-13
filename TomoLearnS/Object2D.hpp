#include <CImg.h>
#include <string>
#include <array>
#include <cstdint>

class Object2D{
public:
	Object2D():objPixSizes{1,1}{  //Default constructor, initialize an empty image
		cimg_image =cimg_library::CImg<uint16_t>();
		cimg_window = cimg_library::CImgDisplay();
	}
	Object2D(std::string imageFilePath, std::array<double, 2> objPixSizes={0.1, 0.1});
	void display(std::string title);
private:
	cimg_library::CImg<uint16_t> cimg_image;
	cimg_library::CImgDisplay cimg_window;
	std::array<double, 2> objPixSizes;
};
