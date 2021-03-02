#pragma once
#include <Eigen/Dense>

enum class FilterType{RamLak,
                      SheppLogan,
                      Cosine,
                      Hamming,
                      Hanning};

class Filter{
public:
	Filter(FilterType filterType, double cutOffIn = 1.0);

	void operator()(Eigen::MatrixXcd& fftOfSinogram);
private:
	FilterType filterType;
	double cutOff;

	void RamLakFilter    (Eigen::MatrixXcd& fftOfSinogram);
	//void SheppLoganFilter(Eigen::MatrixXcd& fftOfSinogram);
	//void CosineFilter    (Eigen::MatrixXcd& fftOfSinogram);
	//void HammingFilter   (Eigen::MatrixXcd& fftOfSinogram);
	//void HanningFilter   (Eigen::MatrixXcd& fftOfSinogram);
};
