#include <TomoLearnS/Gen1CT.hpp>

typedef std::vector<double> Row;
typedef std::vector<Row> Matrix;

Gen1CT::Gen1CT():detWidth{100},pixNum{100},object{nullptr}{
	};

void Gen1CT::putObject(Object2D* sampleObject){
	object = sampleObject;
}

void Gen1CT::measure(const std::vector<double>& angles, int raysPerPixel){
	const int numAngles = angles.size();

	Matrix sinogram(numAngles, Row(pixNum)); //Matrix with (pixNum)x(numAngles) size
}
