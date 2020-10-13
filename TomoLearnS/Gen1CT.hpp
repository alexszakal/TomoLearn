#include <TomoLearnS/Object2D.hpp>
#include <vector>

class Gen1CT{
public:
	Gen1CT();
	Gen1CT(double detWidth, unsigned int pixNum):detWidth{detWidth},pixNum{pixNum} { };
	void putObject(Object2D object);
	void measure(const std::vector<double>& angles, int raysPerPixel=1);
	void filteredBackProjection();

private:
	const double detWidth;
	const unsigned int pixNum;
	Object2D object;
};
