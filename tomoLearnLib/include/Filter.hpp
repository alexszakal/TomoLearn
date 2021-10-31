#pragma once
#include <Eigen/Dense>

enum class FilterType{RamLak,
                      SheppLogan,
                      Cosine,
                      Hann,
                      Hamming};

class Filter{
public:
	Filter(FilterType filterType, double cutOffIn = 1.0);

	void operator()(Eigen::MatrixXcd& fftOfSinogram);
private:
	FilterType filterType;
	double cutOff;

	void RamLakFilter    (Eigen::MatrixXcd& freqFilter);
	void SheppLoganFilter(Eigen::MatrixXcd& freqFilter);
	void CosineFilter    (Eigen::MatrixXcd& freqFilter);
	void HannFilter   (Eigen::MatrixXcd& freqFilter);
	void HammingFilter   (Eigen::MatrixXcd& freqFilter);
};
